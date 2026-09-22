#!/usr/bin/env python3
"""Tabulate, for each degree D in 3..100, the levels N for which X_0(N) is D-elliptic.

This transforms the data in `results/quadratic_forms_x0.csv` into the data in `results/D-elliptic_levels_by_degree.csv`.

Usage: python levels_by_degree.py [input.csv] [output.csv]

The defaults for input.csv output.csv are `results/quadratic_forms_x0.csv` and `results/D-elliptic_levels_by_degree.csv`.
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
                    levels[D].add(int(row["N"]))
    with open(dst, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["degree", "level"])
        for D in DEGREES:
            writer.writerow([D, " ".join(map(str, sorted(levels[D])))])


if __name__ == "__main__":
    main()
