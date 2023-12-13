# tdb: Tandem repeat database and analysis queries
![coverage](imgs/coverage.svg)
![pylint](imgs/pylint.svg)

## Installation

```bash
git clone https://github.com/ACEnglish/tdb.git
cd tdb/
python3 -m pip install . 
```

## Quick Start

Convert a tdb compatible VCF and run a query with:
```
tdb create -o output.tdb input.vcf.gz
tdb query allele_cnts output.tdb
```

## Tutorials

See the [wiki](https://github.com/acenglish/tdb/wiki) or start with
[notebooks/Introduction.ipynb](https://github.com/ACEnglish/tdb/blob/main/notebooks/Introduction.ipynb).
