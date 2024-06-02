#!/usr/bin/env python3
"""
tdb main entrypoint
"""
import argparse

from tdb.create import create_main
from tdb.query import query_main
from tdb.deid import deid_main
from tdb.dump import dump_main
from tdb.merge import merge_main
from tdb.merge_batch import merge_batch_main
from tdb.debug_check import debug_check_main

CMDS = {
    "create": create_main,
    "merge": merge_main,
    "bigmerge": merge_batch_main
    "query": query_main,
    "deid": deid_main,
    "dump": dump_main,
}


cmd_str = ", ".join(list(CMDS.keys()))
USAGE = f"""
tdb v0.2.0 - Tandem repeat database tools

Commands: {cmd_str}"""

def main():
    """
    Main entrypoint for tdb
    """
    # Hidden commands
    CMDS['dbg_eq'] = debug_check_main
    parser = argparse.ArgumentParser(prog="tdb", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=CMDS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    args = parser.parse_args()

    CMDS[args.cmd](args.options)

if __name__ == '__main__':
    main()
