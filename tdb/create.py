"""
Turn a vcf into a tdb
"""
import gc
import os
import sys
import logging
import argparse

import pysam
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import tdb

DTYPES = {"LocusID": pa.uint32(),
          "chrom": pa.string(),
          "start": pa.uint32(),
          "end": pa.uint32(),
          "allele_number": pa.uint16(),
          "allele_length": pa.uint16(),
          "sequence": pa.binary(),
          "spanning_reads": pa.uint16(),
          "length_range_lower": pa.uint16(),
          "length_range_upper": pa.uint16(),
          "average_methylation": pa.float32()}

L_COLUMNS = ["LocusID", "chrom", "start", "end"]
A_COLUMNS = ["LocusID", "allele_number", "allele_length", "sequence"]
S_COLUMNS = ["LocusID", "allele_number", "spanning_reads", "length_range_lower",
             "length_range_upper", "average_methylation"]

L_SCHEMA = pa.schema({key: DTYPES[key] for key in L_COLUMNS})
A_SCHEMA = pa.schema({key: DTYPES[key] for key in A_COLUMNS})
S_SCHEMA = pa.schema({key: DTYPES[key] for key in S_COLUMNS})


def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False

    if os.path.exists(args.output):
        logging.error(f"Output {args.output} already exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"Output {args.output} must end with `.tdb`")
        check_fail = True
    if not os.path.exists(args.input):
        logging.error(f"Input {args.input} does not exist")
        check_fail = True
    if not args.input.rstrip('/').endswith((".vcf", ".vcf.gz")):
        logging.error(
            f"Unrecognized file extension on {args.input}. Expected .vcf .vcf.gz")
        check_fail = True
    return check_fail


def make_parquets(samples, out_dir, compression):
    """
    Parquet writer for the tables
    """
    ret = {}
    comp = "GZIP" if compression else None
    ret['locus'] = pq.ParquetWriter(os.path.join(out_dir, 'locus.pq'),
                                    L_SCHEMA, compression=comp)

    ret['allele'] = pq.ParquetWriter(os.path.join(out_dir, 'allele.pq'),
                                     A_SCHEMA, compression=comp)
    ret['sample'] = {}
    for s in samples:
        fn = os.path.join(out_dir, f"sample.{s}.pq")
        ret['sample'][s] = pq.ParquetWriter(fn, S_SCHEMA, compression=comp)
    return ret


def sample_extract(locus_id, fmt):
    """
    Given a dict from a vcf record sample, turn them into sample rows
    """
    ret = []
    view = zip(fmt['GT'], fmt['SD'], fmt['ALLR'], fmt['AM'])
    for an, sd, allr, am in view:
        if an is None:
            continue
        lrl, lru = map(int, allr.split('-'))
        ret.append([locus_id, an, sd, lrl, lru, am])
    return ret


def translate_entry(entry, locus_id):
    """
    return three things,
    a list of LocusID, chrom, start, end
    a list of allele rows
    a dictionary of sample: list of sample rows
    """
    locus = [locus_id, entry.chrom, entry.start, entry.stop]
    alleles = [(locus_id, allele_number, len(sequence),
                b"" if sequence in [None, "."] else sequence.encode("utf8"))
               for allele_number, sequence in enumerate(entry.alleles)]
    samples = {sample: sample_extract(locus_id, fmt)
               for sample, fmt in entry.samples.items()}
    return locus, alleles, samples


def convert_buffer(vcf, samples, stats, avail_mem):
    """
    Converts a number of vcf entries.
    Tries to monitor memory to not buffer too many
    """
    m_buffer = {'locus': [],
                'allele': [],
                'sample': {s: [] for s in samples}
                }
    # Flag for telling main loop when we're finished
    cvt_any = False
    used_mem = int(avail_mem * 0.20)
    while avail_mem > used_mem:
        try:
            entry = next(vcf)
        except StopIteration:
            break

        cvt_any = True
        cur_l, cur_a, cur_s = translate_entry(entry, stats['locus'])
        m_buffer['locus'].append(cur_l)
        m_buffer['allele'].extend(cur_a)
        num_samples = 0
        for name, rows in cur_s.items():
            num_samples += len(rows)
            m_buffer['sample'][name].extend(rows)

        used_mem += sys.getsizeof(cur_l)
        used_mem += sys.getsizeof(cur_a)
        used_mem += 400 * num_samples

        stats['locus'] += 1
        stats['allele'] += len(cur_a)
        stats['sample'] += num_samples

    return m_buffer, cvt_any


def write_tables(cur_tables, tables):
    """
    Write the cur_tables entries to the output tables
    """
    ldf = pd.DataFrame(cur_tables["locus"], columns=L_COLUMNS, copy=False)
    locus = pa.Table.from_pandas(ldf, schema=L_SCHEMA, preserve_index=False)
    tables['locus'].write(locus)

    adf = pd.DataFrame(cur_tables["allele"], columns=A_COLUMNS, copy=False)
    allele = pa.Table.from_pandas(adf, schema=A_SCHEMA, preserve_index=False)
    tables['allele'].write(allele)

    for name, out_samp in tables["sample"].items():
        sdf = pd.DataFrame(cur_tables["sample"][name],
                           columns=S_COLUMNS, copy=False)
        sample = pa.Table.from_pandas(sdf, schema=S_SCHEMA,
                                      preserve_index=False)
        out_samp.write(sample)


def create_main(args):
    """
    Create a new tdb from multiple input calls
    """
    parser = argparse.ArgumentParser(prog="tdb create", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output tdb directory")
    parser.add_argument("--mem", metavar="MEM", type=int, default=4,
                        help="Memory in GB available to buffer reading (%(default)s)")
    parser.add_argument("--no-compress", action="store_false",
                        help="Don't compress the database")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument("input", metavar="IN",
                        help="VCF file")
    args = parser.parse_args(args)

    tdb.setup_logging(args.debug)

    if check_args(args):
        logging.error("Cannot create database. Exiting")
        sys.exit(1)

    tdb.setup_logging()
    avail_mem = args.mem * 1e9

    os.mkdir(args.output)

    old = pysam.set_verbosity(0)  # suppress non-indexed warning
    vcf = pysam.VariantFile(args.input)
    pysam.set_verbosity(old)  # turn back on
    samples = list(vcf.header.samples)
    stats = {"locus": 0, "allele": 0, "sample": 0}

    tables = make_parquets(samples, args.output, args.no_compress)
    logging.info("Converting VCF with %d samples", len(samples))
    while True:
        cur_tables, cvt_any = convert_buffer(vcf, samples, stats, avail_mem)
        if not cvt_any:
            break
        logging.info("Writing batch. Row totals %s", stats)
        write_tables(cur_tables, tables)
        del cur_tables
        gc.collect()

    tables['locus'].close()
    tables['allele'].close()
    for st in tables['sample'].values():
        st.close()

    logging.info("Finished")
