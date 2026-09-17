#!/usr/bin/env python3
"""Format elliptic-map results, or compute and save them as a text file.

Usage (in the Sage environment): python d_elliptic_output.py 3 [output.txt] [--with-solution]
"""

from sage.all import Integer, PolynomialRing, ZZ, divisors
from pathlib import Path
import argparse
import sys


def format_results(result):
    """Format N, E', L and Q, with a solution only when requested."""
    lines = ["D=%s; status=%s; N <= %s" %
             (result["D"], result["status"], result["N_bound"])]
    for row in result["records"]:
        B, A = row["old_coordinate_basis"], row["target_degree_matrix"]
        rank = A.nrows()
        ring = PolynomialRing(ZZ, names=["z%s" % (i + 1) for i in range(rank)])
        z = ring.gens()
        Q = sum(A[i, j] * z[i] * z[j]
                for i in range(rank) for j in range(rank))
        lines.extend([
            "\nN=%s; M=%s; E=%s; E'=%s" %
            (row["N"], row["M"], row["strong_label"], row["target_label"]),
            "Old-coordinate divisors: %s" % divisors(ZZ(row["R"])),
            "L = B*ZZ^%s, B =" % rank,
            str(B),
            "Q(z) = %s" % Q,
        ])
        if row["solution_z"] is not None:
            lines.extend([
                "z = %s; Q(z) = %s" % (tuple(row["solution_z"]), result["D"]),
                "x = B*z = %s" % (tuple(row["old_coordinates_x"]),),
            ])
    lines.append("\nPositive levels: %s" %
                 sorted({row["N"] for row in result["records"]}))
    if result["errors"]:
        lines.append("Unresolved cases (not negative answers):")
        for error in result["errors"]:
            context = ", ".join("%s=%s" % (key, value) for key, value in error.items()
                                if key not in ("stage", "message"))
            lines.append("  %s: %s: %s" % (context, error["stage"], error["message"]))
    return "\n".join(lines) + "\n"


def save_results(result, filename):
    """Save an already computed result, without overwriting an existing file."""
    text = format_results(result)
    with Path(filename).open("x", encoding="utf-8") as handle:
        handle.write(text)


def main():
    source_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("D", type=int, help="positive map degree")
    parser.add_argument("filename", nargs="?", help="text file (default: results_D<D>.txt)")
    parser.add_argument("--with-solution", action="store_true", help="also find one integral solution")
    args = parser.parse_args()
    if args.D <= 0:
        parser.error("D must be a positive integer")
    filename = Path(args.filename or ("results_D%s.txt" % args.D))
    if filename.exists():
        parser.error("file already exists; choose another filename: %s" % filename)
    if not filename.parent.is_dir():
        parser.error("output directory does not exist: %s" % filename.parent)

    from sage.repl.load import load
    computation = {"_loaded_as_library": True, "Integer": Integer}
    load(str(source_dir / "d_elliptic_search.sage"), computation)
    result = computation["compute_degree"](args.D, with_solution=args.with_solution)
    save_results(result, filename)
    print("Saved: %s (status=%s)" % (filename, result["status"]))
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(int(main()))
