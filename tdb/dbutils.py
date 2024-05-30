"""
Utilities for interacting with a tdb
"""
import os
import sys
import glob
import logging

import pysam
import truvari
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


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


def set_tdb_types(d):
    """
    Sets tdb datatypes of table columns in place
    """
    l_types = {"LocusID": np.uint32,
               "chrom": str,
               "start": np.uint32,
               "end": np.uint32}
    a_types = {"LocusID": np.uint32,
               "allele_number": np.uint16,
               "allele_length": np.uint16,
               "sequence": bytes}

    d['locus'] = d['locus'].astype(l_types)
    d['allele'] = d['allele'].astype(a_types)

    # last I checked, these types couldn't handle nones
    # s_types = {"LocusID": np.uint32,
    #           "allele_number": np.uint16,
    #           "spanning_reads": np.uint16,
    #           "length_range_lower": np.uint16,
    #           "length_range_upper": np.uint16,
    #           "average_methylation": np.float32}

    # for samp, val in d['sample'].items():
    #    d['sample'][samp] = val.astype(s_types)


def write_samples(samples, output):
    """
    Write tdb samples to output folder
    """
    s_schema = pa.schema([('LocusID', pa.uint32()),
                          ('allele_number', pa.uint16()),
                          ('spanning_reads', pa.uint16()),
                          ('length_range_lower', pa.uint16()),
                          ('length_range_upper', pa.uint16()),
                          ('average_methylation', pa.float32())
                          ])
    for sample, value in samples.items():
        o_fn = os.path.join(output, f"sample.{sample}.pq")
        m_table = pa.Table.from_pandas(value)
        n_table = []
        n_names = []
        for col, dtype in zip(s_schema.names, s_schema.types):
            n_table.append(pa.compute.cast(m_table[col], dtype))
            n_names.append(col)
        n_table = pa.Table.from_arrays(n_table, names=n_names)
        writer = pq.ParquetWriter(o_fn, s_schema, compression='gzip')
        writer.write_table(n_table)
        writer.close()


def dump_tdb(data, output):
    """
    Write tdb data to output folder

    WARNING: will overwrite existing data
    """
    if not os.path.exists(output):
        os.mkdir(output)
    pq_fns = get_tdb_filenames(output)
    data['locus'].to_parquet(pq_fns['locus'], index=False, compression='gzip')
    data['allele'].to_parquet(
        pq_fns['allele'], index=False, compression='gzip')
    write_samples(data['sample'], output)


def pull_alleles(data):
    """
    Turn alleles into a table
    """
    alleles = (pd.DataFrame(data["alleles"].to_list(), index=data.index)
               .reset_index()
               .melt(id_vars='hash', value_name="sequence")
               .drop(columns="variable")
               .dropna()
               .set_index('hash'))
    # Full deletions without anchor base(?)
    alleles.loc[alleles["sequence"] == '.', "sequence"] = ""
    alleles["LocusID"] = data["LocusID"]
    alleles["allele_number"] = alleles.groupby(["LocusID"]).cumcount()
    alleles = (alleles.sort_values(["LocusID", "allele_number"])
               .reset_index(drop=True)
               .drop_duplicates(subset=["LocusID", "sequence"]))
    alleles["allele_length"] = alleles["sequence"].str.len()
    alleles["sequence"] = alleles["sequence"].str.encode()
    alleles = (alleles.sort_values(["LocusID", "allele_number"])
               [["LocusID", "allele_number", "allele_length", "sequence"]]
               .reset_index(drop=True))
    return alleles


def pull_saps(data, sample):
    """
    Turn sample allele properties into a table
    """
    # Remove sites with uninformative genotype information
    data = data[~data[f'{sample}_GT'].isin([(None,)])]
    gt = pd.DataFrame(data[f"{sample}_GT"].to_list(), columns=[
                      "GT1", "GT2"], index=data.index)
    span = pd.DataFrame(data[f"{sample}_SD"].to_list(), columns=[
                        "SD1", "SD2"], index=data.index)
    alci = pd.DataFrame(data[f"{sample}_ALLR"].to_list(), columns=[
                        "LR1", "LR2"], index=data.index)
    meth = pd.DataFrame(data[f"{sample}_AM"].to_list(), columns=[
                        "AM1", "AM2"], index=data.index)
    sap = pd.concat([data[["LocusID"]], span, alci, meth, gt], axis=1)

    renamer = {"SD1": "spanning_reads", "SD2": "spanning_reads",
               "LR1": "LR", "LR2": "LR",
               "AM1": "average_methylation", "AM2": "average_methylation",
               "GT1": "allele_number", "GT2": "allele_number"}
    sap = pd.concat([sap[["LocusID", "GT1", "SD1", "LR1", "AM1"]].rename(columns=renamer),
                     sap[["LocusID", "GT2", "SD2", "LR2", "AM2"]].rename(columns=renamer)],
                    axis=0)
    sap[["length_range_lower", "length_range_upper"]
        ] = sap["LR"].str.split('-', expand=True)
    # chrY has None
    sap = sap[~sap["allele_number"].isna()]
    return sap.drop(columns=["LR"]).reset_index(drop=True)


def vcf_to_tdb(vcf_fn):
    """
    Turn a vcf into an in-memory tdb database
    """
    if not os.path.exists(vcf_fn):
        raise RuntimeError(f"input {vcf_fn} does not exist")

    ret = {}
    old = pysam.set_verbosity(0)  # suppressing no-index warning
    data = truvari.vcf_to_df(vcf_fn, with_info=True,
                             with_format=True, alleles=True)
    pysam.set_verbosity(old)

    logging.info("locus count:\t%d", len(data))
    data["LocusID"] = range(len(data))

    ret["locus"] = data[["LocusID", "chrom", "start", "end"]].reset_index(
        drop=True).copy()  # pylint: disable=unsubscriptable-object # pylint/issues/3139

    logging.info("Wrangling alleles")
    allele_df = pull_alleles(data)
    ret["allele"] = allele_df
    logging.info("allele count:\t%d", len(allele_df))

    logging.info("Pulling samples")
    ret["sample"] = {}
    gt_count = 0
    old = pysam.set_verbosity(0)  # suppressing no-index warning
    for sample in pysam.VariantFile(vcf_fn).header.samples:
        ret['sample'][sample] = pull_saps(data, sample)
        gt_count += len(ret['sample'][sample])
    pysam.set_verbosity(old)
    logging.info("genotype count:\t%d", gt_count)
    set_tdb_types(ret)
    return ret
