#!/usr/bin/env python3
"""
Generates coverage for MOOSE and MOOSE-based applications.

Used primarily in coverage generation on CIVET.
"""

import argparse
import os


def parse_args() -> argparse.Namespace:
    """Parse arguments."""
    parser = argparse.ArgumentParser(
        description="Generates coverage for MOOSE and MOOSE-based applications"
    )

    parent = argparse.ArgumentParser(add_help=False)

    action_parser = parser.add_subparsers(dest="action", help="The action to perform.")
    action_parser.required = True

    def add_common_args(parser: argparse.ArgumentParser):
        parser.add_argument(
            "--work-dir",
            type=str,
            default=os.getcwd(),
            help=f"The stateful work directory (default: {os.getcwd()})",
        )

        if (MOOSE_JOBS := os.environ.get("MOOSE_JOBS")) is not None:
            default_jobs = int(MOOSE_JOBS)
            default_jobs_help = f"MOOSE_JOBS={default_jobs}"
        else:
            default_jobs = 4
            default_jobs_help = "4"
        parser.add_argument(
            "--jobs",
            "-j",
            type=int,
            default=default_jobs,
            help=f"Number of parallel jobs (default: {default_jobs_help})",
        )


    def add_search_dirs(parser: argparse.ArgumentParser):
        parser.add_argument("--search-app",
                            action="store_true",
                            help="Add search directories for an app")
        parser.add_argument(
            "--search-dir",
            action="extend",
            nargs=1,
            type=str,
            help="Add a directory to search for coverage."
        )

    initialize_parser = action_parser.add_parser(
        "initialize", parents=[parent], help="Initializes coverage"
    )
    add_common_args(initialize_parser)
    add_search_dirs(initialize_parser)

    return parser.parse_args()

def action_initialize(args):
    print('hello')
    pass


def main():
    args = parse_args()

    action_function = globals()[f"action_{args.action}"]
    action_function(args)

if __name__ == "__main__":
    main()
