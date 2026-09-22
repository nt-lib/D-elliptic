#!/usr/bin/env python3
"""Tabulate, for each degree D in 3..100, the values of M with a saved positive decision.

Usage: python levels_by_degree.py [input.csv] [output.csv]
"""

import csv
import sys

DEGREES = range(3, 101)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "results/quadratic_forms_x0.csv"
    dst = sys.argv[2] if len(sys.argv) > 2 else "results/D-elliptic_levels_by_degree.csv"
    levels = {D: set() for D in DEGREES}
    with open(src, newline="") as f:
        for row in csv.DictReader(f):
            for D in map(int, row["degrees_with_saved_positive_decision"].split()):
                if D in levels:
                    levels[D].add(int(row["M"]))
    with open(dst, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["degree", "level"])
        for D in DEGREES:
            for M in sorted(levels[D]):
                writer.writerow([D, M])


if __name__ == "__main__":
    main()
