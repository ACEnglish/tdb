"""
Faster creation of a tdb

parse each line and write the outputs.
No need to worry about getting crazy with the indexes or joins
"""
import gc
import os
import sys
import logging
import argparse

import pysam
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import tdb

#pylint: disable=global-statement

DTYPES = {"LocusID": (pa.uint32(), np.uint32),
          "chrom": (pa.string(), str),
          "start": (pa.uint32(), np.uint32),
          "end": (pa.uint32(), np.uint32),
          "allele_number": (pa.uint16(), np.uint16),
          "allele_length": (pa.uint16(), np.uint16),
          "sequence": (pa.binary(), bytes),
          "spanning_reads": (pa.uint16(), np.uint16),
          "length_range_lower": (pa.uint16(), np.uint16),
          "length_range_upper": (pa.uint16(), np.uint16),
          "average_methylation": (pa.float32(), np.float32)}

L_COLUMNS = ["LocusID", "chrom", "start", "end"]
A_COLUMNS = ["LocusID", "allele_number", "allele_length", "sequence"]
S_COLUMNS = ["LocusID", "allele_number", "spanning_reads", "length_range_lower",
             "length_range_upper", "average_methylation"]

AVAILMEM = sys.maxsize
# Give 25% overhead since our memory tracking probably underestimates
USEDMEM = int(AVAILMEM * 0.75)

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
        logging.error(f"Unrecognized file extension on {args.input}. Expected .vcf .vcf.gz")
        check_fail = True
    return check_fail

def make_locus_writer(output_file, comp):
    """
    Parquet writer for the locus table
    """
    schema = pa.schema([
        pa.field(key, DTYPES[key][0])  for key in L_COLUMNS
        ])
    return pq.ParquetWriter(output_file, schema, compression=comp)

def make_allele_writer(output_file, comp):
    """
    Parquet writer for the allele table
    """
    schema = pa.schema([
        pa.field(key, DTYPES[key][0]) for key in A_COLUMNS
    ])
    return pq.ParquetWriter(output_file, schema, compression=comp)

def make_sample_writer(output_file, comp):
    """
    Parquet writer for the sample table
    """
    schema = pa.schema([
        pa.field(key, DTYPES[key][0]) for key in S_COLUMNS
    ])
    return pq.ParquetWriter(output_file, schema, compression=comp)

def make_parquets(samples, out_dir, compression):
    """
    Make the parquet file handlers for the output database
    """
    comp = "GZIP" if compression else None
    ret = {}
    ret['locus'] = make_locus_writer(os.path.join(out_dir, 'locus.pq'), comp)
    ret['allele'] = make_allele_writer(os.path.join(out_dir, 'allele.pq'), comp)
    ret['sample'] = {}
    for name in samples:
        ret['sample'][name] = make_sample_writer(os.path.join(out_dir, f"sample.{name}.pq"), comp)
    return ret

def sample_extract(locus_id, fmt_fields):
    """
    Given a dict from a vcf record sample, turn them into sample rows
    """
    ret = []
    view = zip(fmt_fields['GT'],
               fmt_fields['SD'],
               fmt_fields['ALLR'],
               fmt_fields['AM'])
    for an, sd, allr, am in view:
        # None isn't imported
        if an is None:
            continue
        lrl, lru = allr.split('-')
        lrl = int(lrl)
        lru = int(lru)
        ret.append([locus_id, an, sd, lrl, lru, am])
    return ret

def translate_entry(entry, locus_id):
    """
    return three things,
    a list of LocusID, chrom, start, end
    a list of allele rows
    a dictionary of sample: list of sample rows
    """
    global USEDMEM
    locus = [locus_id, entry.chrom, entry.start, entry.stop]
    USEDMEM += sys.getsizeof(locus)
    alleles = [(locus_id, allele_number, len(sequence), sequence.encode("utf8"))
               for allele_number, sequence in enumerate(entry.alleles)]
    USEDMEM += sys.getsizeof(alleles)
    samples = {}
    for sample, m_d in entry.samples.items():
        samples[sample] = sample_extract(locus_id, m_d)
    # Approximate usage of each row of a sample table
    USEDMEM += 400 * len(samples)
    return locus, alleles, samples

def convert_buffer(vcf, samples, stats):
    """
    Converts a number of vcf entries.
    Tries to monitor memory to not buffer too many
    """
    m_buffer = {'locus':[],
              'allele':[],
              'sample': {_:[] for _ in samples}
              }
    # Flag for telling main loop when we're finished
    cvt_any = False
    while AVAILMEM > USEDMEM:
        try:
            entry = next(vcf)
        except StopIteration:
            break

        cvt_any = True
        cur_locus, cur_allele, cur_sample = translate_entry(entry, stats['locus'])

        m_buffer['locus'].append(cur_locus)
        m_buffer['allele'].extend(cur_allele)
        num_samples = 0
        for name, rows in cur_sample.items():
            num_samples += len(rows)
            m_buffer['sample'][name].extend(rows)

        stats['locus'] += 1
        stats['allele'] += len(cur_allele)
        stats['sample'] += num_samples

    return m_buffer, cvt_any

def write_tables(cur_tables, tables):
    """
    Write the cur_tables entries to the output tables
    So since cur_tables will have a list of vales, I should probably do pa.Table from pandas
    """
    global USEDMEM
    schema = pa.schema({key: DTYPES[key][0] for key in L_COLUMNS})
    ldf = pd.DataFrame(cur_tables["locus"], columns=L_COLUMNS, copy=False)
    locus = pa.Table.from_pandas(ldf, schema=schema, preserve_index=False)
    tables['locus'].write(locus)

    schema = pa.schema({key: DTYPES[key][0] for key in A_COLUMNS})
    adf = pd.DataFrame(cur_tables["allele"], columns=A_COLUMNS, copy=False)
    allele = pa.Table.from_pandas(adf, schema=schema, preserve_index=False)
    tables['allele'].write(allele)

    schema = pa.schema({key: DTYPES[key][0] for key in S_COLUMNS})
    for name, out_samp in tables["sample"].items():
        sdf = pd.DataFrame(cur_tables["sample"][name], columns=S_COLUMNS, copy=False)
        sample = pa.Table.from_pandas(sdf, schema=schema, preserve_index=False)
        out_samp.write(sample)
    # Reset memory
    USEDMEM = int(AVAILMEM * 0.75)

def create_main(args):
    """
    Create a new tdb from multiple input calls
    """
    global AVAILMEM
    parser = argparse.ArgumentParser(prog="tdb create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output tdb directory")
    parser.add_argument("--mem", metavar="MEM", type=int,
                        help="Memory in GB available to buffer reading (unlimited)")
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
    if args.mem is not None:
        AVAILMEM = args.mem * 1e9

    os.mkdir(args.output)

    vcf = pysam.VariantFile(args.input)
    samples = list(vcf.header.samples)
    stats = {"locus":0, "allele":0, "sample":0}

    tables = make_parquets(samples, args.output, args.no_compress)
    logging.info("Converting VCF with %d samples", len(samples))
    while True:
        cur_tables, cvt_any = convert_buffer(vcf, samples, stats)
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
