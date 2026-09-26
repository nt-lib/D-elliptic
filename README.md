# Degree-D maps from $X_0(N)$ to elliptic curves

SageMath code accompanying the paper: $D$-elliptic modular curves $X_0(N)$ ([arXiv:2609.24363](https://arxiv.org/abs/2609.24363)). For a positive integer D, it computes all levels N with $genus(X_0(N)) \geq 2$ admitting a degree $D$ map over $Q$ to an
elliptic curve. It also gives the positive definite quadratic degree form whose values give all possible degrees of the map $X_0(N) \to E$. Optionally, it provides an integral solution expressing a degree $D$ map in terms of degeneracy maps.

## Requirements

Use SageMath 10.9 with its Cremona elliptic-curve database. Install `mdsage`
using the Python interpreter in your Sage environment (activate that environment
first if using conda):

```sh
python -m pip install "git+https://github.com/koffie/mdsage.git@602242974e6edc01e7466c46949d40ff589a6fbf"
```

Keep `data/modular_degrees_M9999.csv` in its bundled location. This is input data
(exact Cremona modular degrees through conductor 9999), not a result file.
This small table (about 450 KB) supplies the modular degrees missing from
Sage's mini Cremona database, so the full Cremona database is not required
to look up these values.
Beyond that table, modular degrees use the modular-abelian-variety calculation.
For strong-Weil sources, the code uses c_E=1 when M < 400000, following
[Cremona's optimality and Manin-constant computations](https://johncremona.github.io/ecdata/manin.txt).
For M >= 400000, PARI's `ellweilcurve` computes the Manin constant and checks
optimality using modular symbols; this can be expensive. The shortcut bound
is independent of the installed database size. Missing data or a failed
computation is reported as unresolved, not merely a conductor above the bound.

## Run

```sh
sage d_elliptic_search.sage 3
```

This displays the results without creating a file. To compute and save the
same readable output instead:

```sh
python d_elliptic_output.py 3
```

Add `--with-solution` to either command to find and include one integral solution:

```sh
sage -- d_elliptic_search.sage 3 --with-solution
python d_elliptic_output.py 3 results_D3.txt --with-solution
```

For the `.sage` command, `--` separates Sage's options from the script's arguments.
Run the `.py` command with the Python interpreter in your Sage environment
(activate that environment first if using conda).

The second command creates only `results_D3.txt` in the current directory.
An optional filename can be supplied, for example
`python d_elliptic_output.py 3 my_results.txt`. Existing files are not overwritten.
Both commands use the same formatter in `d_elliptic_output.py`; each command
performs its own computation once. No JSON files or checkpoints are created.
Large degrees may require substantial memory and time; interrupted searches
must be restarted.

## Output

Each positive record gives N, the conductor M, the strong-Weil curve E and
target E' (Cremona labels), a basis matrix B, and the form Q.

- The increasing divisors of N/M index the old coordinates.
- The **columns of B** span L = B Zⁿ.
- Q(z) = zᵀAz is the target degree form; mixed coefficients are 2Aᵢⱼ.
- With `--with-solution`, z satisfies Q(z) = D; x = Bz gives its old coordinates.

For D=3, N=22 and E=E'=`11a1`, one obtains B=I₂,
Q(z)=3z₁²−4z₁z₂+3z₂², and with `--with-solution`, z=(0,1).

The search stops at the first positive target for each strong-Weil pair.
It does not enumerate all targets or solutions, or give rational functions
defining the map. A `partial` result lists unresolved cases; these are not
negative answers. An error-free search reports `complete`.

## Saved results through degree 100

[`x0_elliptic_tables.pdf`](results/x0_elliptic_tables.pdf) presents the levels by degree and positive target degree forms, with Cremona and LMFDB labels, in a readable table.

[`quadratic_forms_x0.csv`](results/quadratic_forms_x0.csv) contains 11,908 distinct positive target-lattice records from completed searches for 3 <= D <= 100. Each row gives N, M, the source and target curve labels (Cremona and LMFDB), the isogeny degree, the differential scaling factor c_u, the lattice basis B, the degree matrix A, its polynomial, and a space-separated list of degrees with a saved positive decision. Matrices are encoded as JSON arrays with exact rational entries. The polynomial variables x1, x2, ... are lattice coordinates (z in the code), not old coordinates Bz.

## Code structure

```text
.
├── README.md
├── LICENSE
├── d_elliptic_search.sage      # computation and console command
├── d_elliptic_output.py        # formatting and file-output command
├── data/
│   ├── modular_degrees_M9999.csv
│   └── LICENSE.ecdata
└── results/
    ├── x0_elliptic_tables.pdf
    └── quadratic_forms_x0.csv
```

Start with `compute_degree` and `evaluate_pair` for the search flow. The numbered
sections in `d_elliptic_search.sage` group each stage with its supporting functions.
The computation uses Sage syntax; the output module uses Python syntax but
still needs Sage for polynomial formatting and for running the computation.

| Stage | Functions and responsibilities |
| --- | --- |
| 1. Ogg sieve | `ogg_levels`: derive the finite bounds and select candidate levels. |
| 2. Sources and filters | `strong_weil_candidates`: enumerate sources; `strong_weil_factor`: construct a factor when needed; `strong_weil_manin_constant`: determine c_E; `isogeny_target_data`: enumerate targets and differential scalars. |
| 2. Divisibility and refined Ogg | `old_coordinate_denominator`: denominator bound; `passes_point_count_filter`: point-count obstruction; `_point_count_level_data` and `_point_count_trace`: cached inputs. |
| 3. Old degree form | `old_degree_form`: compute G using `mdsage.degree_pairing`; `surviving_isogeny_degrees`: old-form obstruction. |
| 4. Lattice selection | `old_coordinate_basis`: select the fixed, Shimura, or homology route and return B. |
| 4. Shimura route | `embedded_shimura_data`: certify applicability; `shimura_factor_data`: cache source data; `embedded_shimura_kernel_equal`: compare embedded kernels; `shimura_kernel_intersection_order`: compute the intersection order. |
| 4. Homology route | `degeneracy_matrices`: degeneracy maps; `integral_preimage_basis`: integral preimage by saturation; `dual_isogeny_lattice`: dual image with real-component correction if needed. |
| 4. Homology support | `integral_homology_star_matrix`: complex conjugation; `modular_symbol_ratio`: signed-symbol ratios; `manin_cusps`: finite generating cusp set. |
| 5. Target form and solution | `evaluate_target`: construct Q and decide representation; `representation_witness`: optionally find one exact integral solution (its internal `search` does the enumeration). |
| 6. Search | `compute_degree`: loop over N and M; `evaluate_pair`: apply the stages to one source (its internal `first_positive` loops over targets); `_error_record`: preserve unresolved failures. |

Shared utilities are `curve_key` and `curve_label` (curve identifiers) and
`load_modular_degrees` (CSV input). `main` handles the command line.
In `d_elliptic_output.py`, `format_results` formats the result,
`save_results` writes it, and `main` runs the file-output command.

## Algorithm

```text
Input: degree D
    |
    v
Finite level bound + Ogg sieve ------------------ reject --> next N
    |
    v
Strong-Weil curves E of conductor M dividing N
    |
    v
Modular-degree divisibility --------------------- reject --> next E
    |
    v
Refined Ogg sieve using #E(F_{p^2}) -------------- reject --> next E
    |
    v
Old degree form G; isogenous targets u: E -> E'
    |
    v
Old-form obstruction --------------------------- reject --> next E'
    |
    v
Surviving E' (fixed-basis targets first)
    |
    +-- c_E * c_u = 1 ---------------------------> B = I
    |                                               |
    +-- c_E * c_u != 1                               |
            |                                       |
            v                                       |
        Cyclic Shimura shortcut certified?          |
        (Ling or direct kernel comparison)          |
            |                                       |
            +-- yes --> Shimura lattice basis B ----+
            |                                       |
            +-- no ---> Integral homology basis B --+
                                                    |
    +-----------------------------------------------+
    |
    v
L = B Z^r; A = deg(u) B^T G B; Q(z) = z^T A z
    |
    v
Does Q represent D? ----------------------------- no --> next E'
    |
   yes
    |
    v
With --with-solution: find z with Q(z) = D
    |
    v
Output N, E', L, Q (optionally z); next E
```

Here c_E is the Manin constant and c_u the differential scaling factor of u.
The first Ogg sieve uses the upper bound D(p+1)^2. Once E is known,
the refined sieve replaces it by D #E(F_{p^2}) = D((p+1)^2 - a_p(E)^2),
checking primes p <= 31 not dividing N. This second Ogg check is applied
before all lattice branches, not only before the homology fallback.
The direct comparison checks the embedded Shimura and degeneracy kernels.
If the old-form filter is unavailable, the code proceeds to the lattice
calculation. If the Shimura criterion does not apply or a required operation
is unavailable, it uses integral homology. Unexpected failures in the Shimura
calculation are reported as unresolved rather than silently triggering that
fallback. Errors in the remaining computation are also reported as unresolved,
not as rejections. After exhausting E and E', the search continues to the
next surviving N.

Lattices and witnesses use exact ZZ/QQ arithmetic. The representation decision
uses PARI `qfminim`, which uses floating-point operations internally; negative
decisions are computational evidence, not independently certified proofs.
With `--with-solution`, positive solutions are obtained by exact enumeration.

The old degree form is computed by `mdsage.degree_pairing`. Homology matrices
use the integral basis of `factor.lattice()`, not the rational modular-symbol
basis. The integral-preimage lattice is computed using SageMath's saturation
operation, which uses HNF internally. The returned basis need not be in
column-HNF form, but spans the same lattice.

## AI assistance

Codex (GPT-5.6) was used to assist in writing and refactoring the code.
The implementation was reviewed against the algorithm in the accompanying
paper, and computational checks were performed.

## License

Code: MIT License (`LICENSE`). Cremona data: Artistic License 2.0
(`data/LICENSE.ecdata`).
