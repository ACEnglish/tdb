#!/usr/bin/env python3
"""
tdb main entrypoint
"""
import argparse

from tdb.create import create_main
from tdb.query import query_main
from tdb.append import append_main
from tdb.deid import deid_main

CMDS = {
    "create": ("Create a tdb", create_main),
    "query": ("Query a tdb", query_main),
    "append": ("Append a VCF/tdb to an existing tdb", append_main),
    "deid": ("Deidentify a tdb", deid_main),
}


cmd_str = "\n".join([f"    {k:9} {t[0]}" for k,t in CMDS.items()])
USAGE = f"""
tdb v0.0.1 - Tandem repeat database tools

Commands:
{cmd_str}"""

def main():
    """
    Main entrypoint for tdb
    """
    parser = argparse.ArgumentParser(prog="tdb", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=CMDS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    args = parser.parse_args()

    CMDS[args.cmd][1](args.options)

if __name__ == '__main__':
    main()
