"""
Utilities for interacting with a tdb
"""
import os
import sys
import glob
import logging
import tempfile
import warnings

import pyarrow.parquet as pq


def setup_logging(debug=False, stream=sys.stderr,
                  log_format="%(asctime)s [%(levelname)s] %(message)s"):
    """
    Create default logger

    :param `debug`: Set log-level to logging.DEBUG
    :type `debug`: boolean, optional
    :param `stream`: Where log is written
    :type `stream`: file handler, optional
    :param `log_format`: Format of log lines
    :type `log_format`: string, optional
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)

    def sendWarningsToLog(message, category, filename, lineno, *args, **kwargs):  # pylint: disable=unused-argument
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno,
                        category.__name__, message)

    warnings.showwarning = sendWarningsToLog


def make_temp_filename(tmpdir=None, suffix=""):
    """
    Get a random filename in a tmpdir with an optional extension
    """
    if tmpdir is None:
        tmpdir = tempfile._get_default_tempdir()  # pylint: disable=protected-access
    fn = os.path.join(tmpdir, next(tempfile._get_candidate_names())) + suffix  # pylint: disable=protected-access
    return fn


def get_tdb_samplenames(file):
    """
    Parses the sample name from a tdb sample parquet files
    """
    ret = []
    full_path = os.path.abspath(os.path.expanduser(file))
    for i in glob.glob(os.path.join(full_path, "sample.*.pq")):
        ret.append(os.path.basename(i)[len('sample.'):-len('.pq')])
    return ret


def get_tdb_filenames(dbname):
    """
    Return names of parquet table files in a tdb
    """
    l_fn = os.path.join(dbname, 'locus.pq')
    a_fn = os.path.join(dbname, 'allele.pq')
    full_path = os.path.abspath(os.path.expanduser(dbname))
    s_files = glob.glob(os.path.join(full_path, 'sample.*.pq'))
    s_names = get_tdb_samplenames(dbname)
    s_dict = dict(zip(s_names, s_files))

    return {'locus': l_fn, 'allele': a_fn, 'sample': s_dict}


def load_tdb(dbname, samples=None, lfilters=None, afilters=None, sfilters=None):
    """
    Loads tdb into DataFrames
    returns dict of {'locus': DataFrame, 'allele': DataFrame, 'sample': {'sname': DataFrame}}

    If samples is provided, only a subset of sample tables are loaded.

    The (l)ocus, (a)llele, and (s)ample filters are passed to pyarrow.parquet.read_table
    filters during loading. Filters allow pulling subsets of data and have structure of
      List[Tuple] or List[List[Tuple]] or None (default)

    From their documentation:
      Each tuple has format: (key, op, value) and compares the key with the value. The
    supported op are: = or ==, !=, <, >, <=, >=, in and not in. If the op is in or not in,
    the value must be a collection such as a list, a set or a tuple.

    If a subset of loci are loaded via lfilters, then a filter of ('LocusID', 'in', loaded_locusids)
    is added to afilters and sfilters

    Example:
        >>> import tdb
        >>> loci = [("LocusID", "in", [23, 32, 77])] # Only load these loci
        >>> alleles = [("allele_number", '>', 0)] # And only their non-ref alleles
        >>> db = tdb.load_tdb("repo_utils/test_files/tdb/merge1.tdb", \
                lfilters=loci, \
                afilters=alleles, \
                samples=["HG02630"])
        >>> len(db['locus'])
        3
        >>> len(db['allele'])
        6
        >>> len(db['sample'])
        1
    """
    def add_filter(filts, n_filt):
        """
        Add a new filter
        """
        if filts is None:
            return n_filt
        if isinstance(filts[0], tuple):
            return n_filt + filts
        filts.insert(0, n_filt)
        return filts
    names = get_tdb_filenames(dbname)
    ret = {}
    ret['locus'] = pq.read_table(names['locus'], filters=lfilters).to_pandas()
    if lfilters:
        loci = [("LocusID", "in", ret['locus']['LocusID'].values)]
        afilters = add_filter(afilters, loci)
        sfilters = add_filter(sfilters, loci)

    ret['allele'] = pq.read_table(
        names['allele'], filters=afilters).to_pandas()

    # backwards compatibility
    if isinstance(ret['allele']['sequence'].iloc[0], str):
        ret['allele']['sequence'] = ret['allele']['sequence'].str.encode(
            'utf-8')

    ret['sample'] = {}
    samp_to_fetch = samples if samples is not None else names["sample"].keys()
    for samp in samp_to_fetch:
        if samp not in names['sample']:
            logging.error(
                "Unable to find sample table `sample.%s.pq` in the tdb", samp)
            sys.exit(1)
        ret['sample'][samp] = pq.read_table(
            names['sample'][samp], filters=sfilters).to_pandas()
    return ret
