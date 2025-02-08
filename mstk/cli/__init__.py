import argparse
from . import analyze_rdf, build, export, ffconv, logplot, topconv, trjconv, wham_pp


def main():
    parser = argparse.ArgumentParser(prog="mstk", description="Molecular Simulation Toolkit",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Add subcommands
    build.add_subcommand(subparsers)
    export.add_subcommand(subparsers)

    ffconv.add_subcommand(subparsers)
    topconv.add_subcommand(subparsers)
    trjconv.add_subcommand(subparsers)

    logplot.add_subcommand(subparsers)
    analyze_rdf.add_subcommand(subparsers)
    wham_pp.add_subcommand(subparsers)

    args = parser.parse_args()
    args.func(args)  # Call the function associated with the subcommand
