import codecs
import os.path
import subprocess
from setuptools import setup

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

VERSION = "0.2.0"

setup(
    name="tdb",
    version=VERSION,
    author="Adam English",
    author_email="",
    url="https://github.com/ACEnglish/tdb",
    packages=["tdb"],
    license="BSD 3-Clause Clear License",
    description="Tandem repeat database and analysis queries",
    include_package_data=True,
    long_description=open("README.md", encoding="UTF-8").read(),
    long_description_content_type="text/markdown",
    entry_points={
      "console_scripts": [
         "tdb = tdb.__main__:main"
      ]
    },
    install_requires=[
        "duckdb>=0.6.1",
        "numpy>=1.26.4",
        "pandas>=2.2.2",
        "pysam>=0.22.1",
        "pyarrow>=16.1",
    ],
)
