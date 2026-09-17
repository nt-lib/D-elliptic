#!/usr/bin/env sage
"""Search for degree-D maps X_0(N) -> E' over QQ, with genus(X_0(N)) >= 2."""

from sage.all import (
    CremonaDatabase, EllipticCurve, Gamma0, Infinity, Integer, J0, QQ,
    QuadraticForm, ZZ, cached_function, ceil, divisors,
    factor, floor, gcd, identity_matrix, matrix, next_prime,
    number_of_divisors, pari, prime_range, vector,
)
import argparse
import csv
from collections import namedtuple
from pathlib import Path
from types import SimpleNamespace
import sys
import warnings

# Sage's launcher can leave an optional argument separator in sys.argv.
if sys.argv[0] == "--":
    del sys.argv[0]
SOURCE_DIR = Path(sys.argv[0]).resolve().parent

from mdsage import degree_pairing, modular_symbol_elliptic_curves, product_isogeny_map

CREMONA_DATABASE = CremonaDatabase()
CREMONA_BOUND = ZZ(CREMONA_DATABASE.largest_conductor())
# Conservative range with established optimal-curve identification; exclusive.
MANIN_CONSTANT_BOUND = ZZ(400000)

StrongWeilCandidate = namedtuple(
    "StrongWeilCandidate", "source_key label factor curve")
IsogenyTarget = namedtuple(
    "IsogenyTarget", "target_key label curve degree c_u c_dual")


def load_modular_degrees():
    """Read the bundled Cremona modular degrees through conductor 9999."""
    values = {}
    with (SOURCE_DIR / "data" / "modular_degrees_M9999.csv").open() as handle:
        for row in csv.DictReader(handle):
            label, degree = row["cremona_label"], ZZ(row["modular_degree"])
            if not label or degree <= 0 or label in values:
                raise ArithmeticError("invalid modular-degree table")
            values[label] = degree
    if len(values) != 38042:
        raise ArithmeticError("incomplete modular-degree table")
    return values


EXACT_MODULAR_DEGREES = load_modular_degrees()


POINT_COUNT_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


# Step 1. Ogg sieve -----------------------------------------------------------

def ogg_levels(D):
    """Return exhaustive Ogg survivors and the automatically chosen bounds."""
    D = ZZ(D)
    if D <= 0:
        raise ValueError("D must be a positive integer")
    primorial = ZZ(1)
    possible_bounds = []
    p = ZZ(2)
    while True:
        bound = floor(12 * (D * (p + 1)^2 - 1) / (p - 1))
        if primorial > bound:
            N_bound, first_impossible_prime = max(possible_bounds), p
            break
        possible_bounds.append(ZZ(bound))
        primorial *= p
        p = next_prime(p)
    eligible = []
    max_index = ZZ(0)
    for N in range(1, int(N_bound) + 1):
        group = Gamma0(N)
        genus = ZZ(group.genus())
        if genus < 2:
            continue
        index = ZZ(group.index())
        two_to_omega = 2^len(ZZ(N).prime_divisors())
        eligible.append((ZZ(N), genus, index, two_to_omega))
        max_index = max(max_index, index)

    # Larger primes cannot impose a stronger Ogg bound.
    threshold = max(ZZ(11), ceil(QQ(max_index) / (6 * D)))
    prime_bound = ZZ(next_prime(threshold))
    primes = tuple(ZZ(p) for p in prime_range(2, prime_bound + 1))

    survivors = []
    for N, genus, index, two_to_omega in eligible:
        if all(
            N % p == 0 or
            D * (p + 1)^2 >= (p - 1) * index / 12 + two_to_omega
            for p in primes
        ):
            survivors.append((N, genus))
    return {
        "eligible_count": len(eligible),
        "survivors": survivors,
        "N_bound": int(N_bound),
        "prime_bound": int(prime_bound),
        "first_impossible_least_prime": int(first_impossible_prime),
    }


# Step 2. Strong-Weil curves and filters --------------------------------------

@cached_function
def strong_weil_candidates(M):
    """Return keyed strong-Weil curves, delaying database factor creation."""
    M = ZZ(M)
    if M <= CREMONA_BOUND:
        factor_curves = (
            (None, curve)
            for curve in CREMONA_DATABASE.iter_optimal([M])
        )
    else:
        factor_curves = (
            (factor, factor.elliptic_curve())
            for factor in modular_symbol_elliptic_curves(M)
        )
    rows = []
    for factor, curve in factor_curves:
        if ZZ(curve.conductor()) != M:
            raise ArithmeticError("invalid strong-Weil curve at conductor %s" % M)
        if factor is not None and (
                factor.dimension() != 1 or ZZ(factor.level()) != M):
            raise ArithmeticError("invalid strong-Weil factor at conductor %s" % M)
        rows.append(StrongWeilCandidate(
            curve_key(curve), curve_label(curve), factor, curve))
    rows.sort(key=lambda candidate: (candidate.label, candidate.source_key))
    if len({candidate.source_key for candidate in rows}) != len(rows):
        raise ArithmeticError("duplicate strong-Weil factor at conductor %s" % M)
    return tuple(rows)


@cached_function
def strong_weil_factor(source_key):
    """Materialize one strong-Weil factor, only after it is actually needed."""
    E = EllipticCurve(list(source_key))
    M = ZZ(E.conductor())
    if M <= CREMONA_BOUND:
        factor = E.abelian_variety()
    else:
        matches = [candidate.factor for candidate in strong_weil_candidates(M)
                   if candidate.source_key == source_key]
        if len(matches) != 1 or matches[0] is None:
            raise ArithmeticError("the strong-Weil factor was not found uniquely")
        factor = matches[0]
    if factor.dimension() != 1 or ZZ(factor.level()) != M:
        raise ArithmeticError("invalid strong-Weil factor at conductor %s" % M)
    return factor


@cached_function
def strong_weil_manin_constant(source_key):
    """Use c_E=1 below 400000; otherwise compute it by modular symbols."""
    curve = EllipticCurve(list(source_key))
    M = ZZ(curve.conductor())
    if M < MANIN_CONSTANT_BOUND:
        if M <= CREMONA_BOUND and not curve.is_isomorphic(curve.optimal_curve()):
            raise ArithmeticError("the candidate is not the strong-Weil curve")
        return ZZ(1)

    curves, smith_invariants = pari.ellweilcurve(pari.ellinit(list(curve.ainvs())))
    for model, invariants in zip(curves, smith_invariants):
        if not curve.is_isomorphic(EllipticCurve([QQ(a) for a in list(model)[:5]])):
            continue
        if len(invariants) != 2 or invariants[0] != invariants[1]:
            raise ArithmeticError("modular symbols identify a different optimal curve")
        c_E = ZZ(invariants[0])
        if c_E <= 0:
            raise ArithmeticError("the computed Manin constant is not positive")
        return c_E
    raise ArithmeticError("the source is missing from PARI's isogeny class")


def old_coordinate_denominator(c, R):
    r"""Return D_c(R): 1 if R=1, otherwise prod l^(v_l(c)+floor(v_l(R)/2))."""
    c, R = ZZ(c), ZZ(R)
    if c <= 0 or R <= 0:
        raise ValueError("c and R must be positive")
    if R == 1 or c == 1:
        return ZZ(1)
    value = ZZ(1)
    for ell, exponent in factor(c):
        value *= ZZ(ell)^(ZZ(exponent) + R.valuation(ell) // 2)
    return value


def passes_point_count_filter(D, N, source_curve):
    """Test a necessary degree-D condition for a curve of conductor dividing N."""
    D, N = ZZ(D), ZZ(N)
    index, cusp_bound = _point_count_level_data(N)
    source_key = curve_key(source_curve)
    for p in POINT_COUNT_PRIMES:
        if N % p == 0:
            continue
        ap = _point_count_trace(source_key, p)
        if (p - 1) * index + 12 * cusp_bound > 12 * D * ((p + 1)^2 - ap^2):
            return False
    return True


@cached_function
def _point_count_level_data(N):
    return ZZ(Gamma0(N).index()), 2^len(N.prime_divisors())


@cached_function
def _point_count_trace(source_key, p):
    return ZZ(EllipticCurve(list(source_key)).ap(p))


def curve_key(E):
    """Return a database-independent key for an integral elliptic curve model."""
    return tuple(ZZ(a) for a in E.ainvs())


def curve_label(E):
    """Return a readable curve label, falling back to its integral a-invariants."""
    try:
        return str(E.cremona_label())
    except (AttributeError, LookupError, RuntimeError, ValueError):
        return "ainvs:" + ",".join(str(a) for a in curve_key(E))


@cached_function
def isogeny_target_data(source_key):
    """Return all targets and differential data for an integral minimal source."""
    source_curve = EllipticCurve(list(source_key))
    isogeny_class = source_curve.isogeny_class()
    curves = list(isogeny_class.curves)
    degrees = isogeny_class.matrix()
    source = next(
        (i for i, curve in enumerate(curves)
         if tuple(curve.ainvs()) == tuple(source_curve.ainvs())),
        None,
    )
    if source is None:
        raise ValueError("the source must be an integral minimal model")
    # Compose differential scalars along a minimal prime-isogeny path.
    prime_isogenies = isogeny_class.isogenies()
    rows = []
    for j, target in enumerate(curves):
        degree = ZZ(degrees[source, j])
        if degree <= 0:
            raise ArithmeticError("an isogeny-class degree must be positive")
        current, c_u = source, ZZ(1)
        while current != j:
            step = next(
                ((k, morphism) for k, morphism in
                 enumerate(prime_isogenies[current])
                 if morphism and
                 ZZ(morphism.degree()) * degrees[k, j] == degrees[current, j]),
                None)
            if step is None:
                raise ArithmeticError("the prime-isogeny graph is incomplete")
            current, morphism = step
            c_u *= abs(ZZ(morphism.scaling_factor()))
        c_dual, remainder = degree.quo_rem(c_u)
        if remainder:
            raise ArithmeticError(
                "the differential factor does not divide the isogeny degree")
        rows.append(IsogenyTarget(
            curve_key(target), curve_label(target), target,
            degree, c_u, c_dual))
    rows.sort(key=lambda target: (
        target.label, target.target_key, target.degree))
    if not rows or not any(
            target.curve.is_isomorphic(source_curve) and target.degree == 1
            for target in rows):
        raise ArithmeticError("isogeny class is missing its identity target")
    return tuple(rows)


# Step 3. Old degree form and obstruction -------------------------------------

def old_degree_form(source_key, N, modular_degree):
    """Return the old degree form computed by mdsage."""
    source_curve = EllipticCurve(list(source_key))
    N, modular_degree = ZZ(N), ZZ(modular_degree)
    M = ZZ(source_curve.conductor())
    if N <= 0 or N % M:
        raise ValueError("N must be a positive multiple of the conductor")
    if modular_degree <= 0:
        raise ValueError("the modular degree must be positive")

    # Supply the known degree without changing the cached Sage factor.
    modular_factor = strong_weil_factor(source_key)
    factor_view = SimpleNamespace(
        dimension=modular_factor.dimension, level=modular_factor.level,
        modular_symbols=modular_factor.modular_symbols,
        modular_degree=lambda: modular_degree)
    mdsage_form = matrix(ZZ, degree_pairing(factor_view, N))
    if (not mdsage_form.is_symmetric() or
            not mdsage_form.is_positive_definite()):
        raise ArithmeticError("the old degree form is not positive definite symmetric")
    return mdsage_form


def surviving_isogeny_degrees(D, denominator_factor, old_form, target_rows):
    """Return the isogeny degrees surviving the old-form obstruction."""
    D, denominator_factor = ZZ(D), ZZ(denominator_factor)
    if D <= 0 or denominator_factor <= 0:
        raise ValueError("D and the denominator factor must be positive")
    degrees = sorted({ZZ(target.degree) for target in target_rows})
    if not degrees:
        return set()
    values = {
        degree: degree * denominator_factor^2 * D
        for degree in degrees
    }
    bound = max(values.values())
    quadratic_form = QuadraticForm(ZZ, 2 * matrix(ZZ, old_form))
    representation_counts = quadratic_form.representation_number_list(
        int(bound) + 1)
    return {degree for degree in degrees
            if representation_counts[int(values[degree])]}


# Step 4. Old-coordinate lattice ----------------------------------------------

def old_coordinate_basis(
        source_curve, target_curve, degree, c_E, c_u, c_dual, rank,
        shimura=None, homology_matrices=None):
    r"""Return a column basis of the old-coordinate lattice L_{E',N}."""
    degree = ZZ(degree)
    c_E, c_u, c_dual = ZZ(c_E), ZZ(c_u), ZZ(c_dual)
    if min(degree, c_E, c_u, c_dual) <= 0:
        raise ValueError("isogeny and differential data must be positive")
    if c_u * c_dual != degree:
        raise ArithmeticError("c_u*c_dual is not the isogeny degree")

    rank = int(rank)
    if rank <= 0:
        raise ValueError("the old-coordinate rank must be positive")
    if c_E * c_u == 1:
        basis = identity_matrix(QQ, rank)
        method = "fixed"
    elif shimura is not None and shimura["applies"]:
        intersection_order = shimura_kernel_intersection_order(
            shimura, curve_key(source_curve), curve_key(target_curve),
            degree, c_dual)
        basis = matrix(QQ, rank, rank)
        basis[0, 0] = 1
        for j in range(1, rank):
            basis[0, j] = 1 / intersection_order
            basis[j, j] = -1 / intersection_order
        method = "shimura"
    else:
        if homology_matrices is None:
            raise ArithmeticError("missing degeneracy homology matrices")
        if len(homology_matrices) != rank:
            raise ArithmeticError("wrong number of degeneracy homology matrices")
        image = dual_isogeny_lattice(
            curve_key(source_curve), curve_key(target_curve), degree, c_dual)
        image_inverse = matrix(QQ, image).inverse()
        columns = []
        for morphism_matrix in homology_matrices:
            lift = QQ(degree) * matrix(QQ, morphism_matrix) * image_inverse
            columns.append(vector(QQ, lift.list()))
        coefficient_map = matrix(QQ, columns).transpose()
        basis = integral_preimage_basis(coefficient_map)
        method = "hnf"
    basis = matrix(QQ, basis)
    return basis, method


def embedded_shimura_data(source_key, N):
    """Return cyclic Shimura data certified by theorem or direct comparison."""
    N = ZZ(N)
    M = ZZ(EllipticCurve(list(source_key)).conductor())
    if N == M:
        return {
            "applies": True,
            "reason": "N equals M",
            "order": ZZ(1),
            "generator": None,
            "shimura_invariants": (),
        }
    R = N // M
    ling_applies = (R.is_squarefree() and gcd(M, R) == 1 and
                    (M % 2 == 1 or R.is_prime()))
    cached = shimura_factor_data(source_key)
    if not cached["applies"]:
        return cached
    if ling_applies:
        # Ling, Theorem 4(i): intersect Sigma(M)_0 with E^tau(R).
        reason = "certified by Ling's theorem"
    else:
        equal, kernel_invariants = embedded_shimura_kernel_equal(
            cached["factor"], cached["subgroup"], N)
        if not equal:
            return {
                "applies": False,
                "reason": ("direct embedded Shimura-kernel comparison "
                           "is unequal"),
                "shimura_invariants": cached["shimura_invariants"],
                "kernel_invariants": kernel_invariants,
            }
        reason = "verified by direct embedded Shimura-kernel comparison"
    result = dict(cached)
    result["reason"] = reason
    return result


@cached_function
def shimura_factor_data(source_key):
    """Cache the N-independent cyclic Shimura data of a strong Weil curve."""
    factor = strong_weil_factor(source_key)
    subgroup = factor.shimura_subgroup()
    invariants = tuple(ZZ(n) for n in subgroup.invariants())
    if len(invariants) > 1:
        return {"applies": False, "reason": "Shimura subgroup is not cyclic"}
    order = invariants[0] if invariants else ZZ(1)
    generators = [g for g in subgroup.gens()
                  if not g.is_zero() and ZZ(g.additive_order()) == order]
    ambient_generator = (vector(QQ, generators[0].element())
                         if order > 1 and generators else None)
    if order > 1 and ambient_generator is None:
        return {"applies": False, "reason": "no cyclic Shimura generator"}
    generator = None
    if order > 1:
        factor_basis = matrix(QQ, factor.lattice().basis_matrix())
        generator = factor_basis.transpose().solve_right(ambient_generator)
    return {
        "applies": True,
        "factor": factor,
        "subgroup": subgroup,
        "order": order,
        "generator": generator,
        "shimura_invariants": invariants,
    }


def embedded_shimura_kernel_equal(factor, subgroup, N):
    """Compare the actual embedded subgroups, not just their invariant factors."""
    N = ZZ(N)
    R = N // ZZ(factor.level())
    rank = int(number_of_divisors(R))
    xi = product_isogeny_map(factor, N)
    ambient, kernel = xi.domain(), xi.kernel()[0]
    ambient_degree = int(ambient.lattice().degree())
    if ambient_degree % rank:
        raise ArithmeticError("invalid product-isogeny coordinate blocks")
    block = ambient_degree // rank

    generators = []
    for position in range(1, rank):
        for generator in subgroup.gens():
            coordinates = vector(QQ, generator.element())
            if len(coordinates) != block:
                raise ArithmeticError(
                    "unexpected embedded Shimura coordinate length")
            embedded = vector(QQ, ambient_degree)
            embedded[0:block] = -coordinates
            embedded[position * block:(position + 1) * block] = coordinates
            generators.append(embedded)
    expected = (ambient.finite_subgroup(
        generators, field_of_definition=subgroup.field_of_definition())
        if generators else ambient.zero_subgroup())
    if not expected.is_subgroup(kernel):
        raise ArithmeticError(
            "the embedded Shimura subgroup is not in the degeneracy kernel")
    return kernel.is_subgroup(expected), tuple(
        ZZ(value) for value in kernel.invariants())


def shimura_kernel_intersection_order(
        shimura, source_key, target_key, degree, c_dual):
    """Return ``#(ker(u) intersect Sigma_E)`` for a cyclic Shimura group."""
    shimura_order = ZZ(shimura["order"])
    if shimura_order == 1:
        return ZZ(1)
    dual_image = dual_isogeny_lattice(
        source_key, target_key, degree, c_dual).row_module()
    generator = shimura["generator"]
    intersection_order = ZZ(sum(
        ZZ(degree) * k * generator in dual_image
        for k in range(int(shimura_order))))
    if (intersection_order <= 0 or
            shimura_order % intersection_order or
            ZZ(degree) % intersection_order):
        raise ArithmeticError("invalid Shimura-kernel intersection order")
    return intersection_order


def degeneracy_matrices(source_key, N):
    """Return homology matrices shared by one pair's targets, not globally cached."""
    modular_factor = strong_weil_factor(source_key)
    M = ZZ(modular_factor.level())
    JN, JM = J0(ZZ(N)), J0(M)
    projection = JM.projection(modular_factor)
    old_divisors = tuple(sorted(divisors(ZZ(N) // M)))
    matrices = tuple((projection * JN.degeneracy_map(M, d)).matrix()
                     for d in old_divisors)
    return old_divisors, matrices


def integral_preimage_basis(linear_map):
    """Return a column basis of {x in QQ^n : linear_map*x in ZZ^m}, by saturation."""
    linear_map = matrix(QQ, linear_map)
    if linear_map.rank() != linear_map.ncols():
        raise ArithmeticError("coefficient-to-homology map is not injective")
    denominator = linear_map.denominator()
    integral_map = matrix(ZZ, denominator * linear_map)
    saturated = integral_map.column_module().saturation()
    saturated_basis = saturated.basis_matrix().transpose()
    basis = integral_map.change_ring(QQ).solve_right(
        QQ(denominator) * saturated_basis.change_ring(QQ))
    return basis


@cached_function
def dual_isogeny_lattice(source_key, target_key, degree, dual_scaling):
    """Return the dual image in integral homology (row-HNF basis)."""
    degree, scaling = ZZ(degree), ZZ(dual_scaling)
    star = integral_homology_star_matrix(source_key)
    plus_index = scaling * modular_symbol_ratio(
        source_key, target_key, 1)
    minus_index = scaling * modular_symbol_ratio(
        source_key, target_key, -1)
    eigenspaces = {
        sign: ((star.transpose() - sign * identity_matrix(ZZ, 2))
               .right_kernel().intersection(ZZ^2))
        for sign in (1, -1)
    }
    if any(module.rank() != 1 for module in eigenspaces.values()):
        raise ArithmeticError("star eigenspace is not one-dimensional")

    original_error = None
    for corrected in (False, True):
        ratio = QQ(1)
        if corrected:
            source = EllipticCurve(list(source_key))
            target = EllipticCurve(list(target_key))
            ratio = QQ(target.real_components()) / source.real_components()
            if ratio == 1:
                raise original_error
        plus, minus = plus_index * ratio, minus_index * ratio
        try:
            if (plus <= 0 or minus <= 0 or
                    plus.denominator() != 1 or
                    minus.denominator() != 1):
                raise ArithmeticError("nonintegral dual-isogeny eigenspace indices")
            plus, minus = ZZ(plus), ZZ(minus)

            matches = []
            for a in divisors(degree):
                c = degree // a
                for b in range(int(c)):
                    lattice = matrix(ZZ, [[a, b], [0, c]])
                    lattice_module = lattice.row_module()
                    is_stable = all(row * star in lattice_module
                                    for row in lattice_module.basis())
                    if (is_stable and
                            lattice_module.intersection(eigenspaces[1]).index_in(
                                eigenspaces[1]) == plus and
                            lattice_module.intersection(eigenspaces[-1]).index_in(
                                eigenspaces[-1]) == minus):
                        matches.append(lattice)
            if len(matches) != 1:
                raise ArithmeticError("dual-isogeny lattice is not uniquely determined")
            return matches[0]
        except ArithmeticError as error:
            if corrected:
                raise ArithmeticError("%s; real-component correction: %s" %
                                      (original_error, error)) from error
            original_error = error


@cached_function
def integral_homology_star_matrix(source_key):
    """Return star on factor.lattice(), not the rational modular-symbol basis."""
    factor = strong_weil_factor(source_key)
    ambient = factor.modular_symbols().ambient_module()
    cuspidal_integral = ambient.cuspidal_submodule().integral_structure()
    ambient_star_QQ = matrix(QQ,
        ambient.star_involution().matrix().restrict(cuspidal_integral))
    if ambient_star_QQ.denominator() != 1:
        raise ArithmeticError("star does not preserve ambient integral homology")
    ambient_star = matrix(ZZ, ambient_star_QQ)
    star_QQ = matrix(QQ, ambient_star.restrict(factor.lattice()))
    if star_QQ.denominator() != 1:
        raise ArithmeticError("star does not preserve factor integral homology")
    star = matrix(ZZ, star_QQ)
    if star.dimensions() != (2, 2):
        raise ArithmeticError("invalid star involution on elliptic homology")
    return star


def modular_symbol_ratio(source_key, target_key, sign):
    """Return the exact signed-symbol ratio, using all finite Manin endpoints."""
    source_curve = EllipticCurve(list(source_key))
    target_curve = EllipticCurve(list(target_key))
    source_symbol = source_curve.modular_symbol(
        sign=sign, implementation="eclib")
    target_symbol = target_curve.modular_symbol(
        sign=sign, implementation="eclib")
    for cusp in manin_cusps(source_key):
        target_value = target_symbol(cusp)
        if target_value:
            return QQ(source_symbol(cusp)) / QQ(target_value)
    raise ArithmeticError("the sign-%s target modular symbol vanished on "
                          "all Manin generators" % sign)


@cached_function
def manin_cusps(source_key):
    """Return finite cusps coming from a Manin-symbol generating set."""
    E = EllipticCurve(list(source_key))
    ambient = E.modular_symbol_space(sign=0).ambient_module()
    cusps = {
        QQ(endpoint)
        for symbol in ambient.manin_symbols().manin_symbol_list()
        for endpoint in symbol.endpoints()
        if endpoint != Infinity
    }
    if not cusps:
        raise ArithmeticError("the Manin-symbol generating set has no finite cusp")
    return tuple(sorted(cusps))


# Step 5. Target degree form and solution -------------------------------------

def evaluate_target(D, N, candidate, target, modular_degree, c_E, old_form,
                    *, shimura=None, homology_matrices=None, with_solution=False):
    """Build the target degree form and decide whether it represents D."""
    source_curve, target_curve = candidate.curve, target.curve
    D, N, degree = ZZ(D), ZZ(N), ZZ(target.degree)
    modular_degree = ZZ(modular_degree)
    conductor = ZZ(source_curve.conductor())
    if D <= 0 or modular_degree <= 0:
        raise ValueError("D and the modular degree must be positive")
    if N <= 0 or N % conductor:
        raise ValueError("N must be a positive multiple of the conductor")
    c_E, c_u, c_dual = ZZ(c_E), ZZ(target.c_u), ZZ(target.c_dual)
    # 4. Construct the lattice for this target.
    basis, method = old_coordinate_basis(
        source_curve, target_curve, degree, c_E, c_u, c_dual,
        old_form.nrows(),
        shimura=shimura, homology_matrices=homology_matrices)

    # 5. Decide representation; optionally find a degree-D solution.
    target_form_rational = (
        QQ(degree) * basis.transpose() * old_form.change_ring(QQ) * basis)
    if target_form_rational.denominator() != 1:
        raise ArithmeticError("target degree form is not integral")
    target_form = matrix(ZZ, target_form_rational)
    if not target_form.is_symmetric() or not target_form.is_positive_definite():
        raise ArithmeticError("target degree form is not positive definite symmetric")
    # PARI evaluates x^T*target_form*x directly.  Storing zero vectors keeps
    # this single-value test independent of the number of short vectors.
    _, maximum_norm, _ = pari.qfminim(target_form, D, 0)
    represents_D = ZZ(maximum_norm) == D
    solution = None
    old_coordinates = None
    if represents_D and with_solution:
        solution = representation_witness(target_form, D)
        if solution is None:
            raise ArithmeticError("positive representation test has no exact witness")
        x = basis * vector(ZZ, solution)
        old_coordinates = list(x)
    return {
        "D": int(D),
        "N": int(N),
        "M": int(conductor),
        "R": int(N // conductor),
        "strong_label": curve_label(source_curve),
        "target_label": curve_label(target_curve),
        "delta": int(degree),
        "modular_degree": int(modular_degree),
        "c_E": int(c_E),
        "c_u": int(c_u),
        "c_dual": int(c_dual),
        "method": method,
        "represents_D": represents_D,
        "old_coordinate_basis": basis,
        "target_degree_matrix": target_form,
        "old_degree_matrix": old_form,
        "solution_z": solution,
        "old_coordinates_x": old_coordinates,
    }


def representation_witness(form, degree):
    """Find one solution by exact LLL/LDL enumeration; form must be positive definite."""
    A = matrix(ZZ, form)
    D = ZZ(degree)
    if D <= 0 or not A.is_square() or A.nrows() == 0:
        raise ValueError("a nonempty positive-definite form and positive degree are required")
    U = matrix(ZZ, A.LLL_gram())
    reduced = U.transpose() * A * U
    rank = A.nrows()
    _, lower, diagonal_matrix = reduced.change_ring(QQ).block_ldlt(classical=True)
    upper = lower.transpose()
    diagonal = diagonal_matrix.diagonal()
    coordinates = [ZZ(0)] * rank

    def search(i, remaining):
        if i < 0:
            return vector(ZZ, coordinates) if remaining == 0 else None
        shift = sum(upper[i, j] * coordinates[j] for j in range(i + 1, rank))
        shift = QQ(shift)
        p, q = shift.numerator(), shift.denominator()
        radius = ZZ((remaining * q^2 / diagonal[i]).floor()).isqrt()
        lower = ZZ((QQ(-radius - p) / q).ceil())
        upper_bound = ZZ((QQ(radius - p) / q).floor())
        for value in range(int(lower), int(upper_bound) + 1):
            coordinates[i] = ZZ(value)
            rest = remaining - diagonal[i] * (value + shift)^2
            if rest < 0:
                raise ArithmeticError("invalid exact enumeration bound")
            found = search(i - 1, rest)
            if found is not None:
                return found
        return None

    found = search(rank - 1, QQ(D))
    if found is None:
        return None
    z = U * found
    if next(entry for entry in z if entry) < 0:
        z = -z
    return list(z)


# Step 6. Search --------------------------------------------------------------

def compute_degree(D, *, with_solution=False):
    """Return positive records and errors from the finite degree-D search."""
    try:
        degree = ZZ(D)
        if isinstance(D, bool) or QQ(D) != degree or degree <= 0:
            raise ValueError
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("D must be a positive integer") from error
    # 1. Select levels by the Ogg sieve.
    sieve = ogg_levels(degree)
    records, errors = [], []
    for N, genus in sieve["survivors"]:
        # 2. Enumerate sources; evaluate_pair applies stages 2-5 to each.
        for M in divisors(N):
            try:
                candidates = strong_weil_candidates(M)
            except Exception as error:
                errors.append(_error_record("strong-Weil enumeration", error, N=N, M=M))
                continue
            for candidate in candidates:
                positive, unresolved = evaluate_pair(
                    degree, N, M, candidate, with_solution=with_solution)
                records.extend(positive)
                errors.extend(unresolved)
    return {
        "D": degree, "N_bound": sieve["N_bound"],
        "status": "partial" if errors else "complete",
        "records": records, "errors": errors,
    }


def evaluate_pair(D, N, M, candidate, *, with_solution=False):
    """Return positive targets and unresolved errors for one strong-Weil pair."""
    D, N, M = ZZ(D), ZZ(N), ZZ(M)
    label, source_curve = candidate.label, candidate.curve
    source_key = candidate.source_key
    errors = []
    # 2. Determine c_E, then apply divisibility and refined Ogg filters.
    try:
        c_E = strong_weil_manin_constant(source_key)
    except Exception as error:
        return [], [_error_record("Manin constant", error, N=N, M=M, E=label)]
    try:
        modular_degree = EXACT_MODULAR_DEGREES.get(label)
        if modular_degree is None:
            modular_factor = candidate.factor
            if modular_factor is None:
                modular_factor = strong_weil_factor(source_key)
            modular_degree = ZZ(modular_factor.modular_degree())
        denominator_factor = old_coordinate_denominator(c_E, N // M)
        if (denominator_factor * D) % modular_degree:
            return [], errors
        if not passes_point_count_filter(D, N, source_curve):
            return [], errors
        # 3. Build the old degree form and apply its obstruction.
        target_data = isogeny_target_data(source_key)
        old_form = old_degree_form(source_key, N, modular_degree)
    except Exception as error:
        return [], [_error_record("pair setup", error, N=N, M=M, E=label)]

    try:
        represented_degrees = surviving_isogeny_degrees(
            D, denominator_factor, old_form, target_data)
        surviving_targets = [target for target in target_data
                             if target.degree in represented_degrees]
    except Exception as error:
        # The old-form obstruction is optional; the lattice calculation remains.
        warnings.warn("N=%s, E=%s: old-form filter bypassed (%s)" % (N, label, error))
        surviving_targets = list(target_data)

    # 4-5. Try the fixed lattice first; prepare other lattices only if needed.
    fixed_targets = [target for target in surviving_targets if c_E * target.c_u == 1]
    lattice_targets = [target for target in surviving_targets if c_E * target.c_u != 1]

    def first_positive(rows, shimura=None, homology_matrices=None):
        for target in rows:
            try:
                result = evaluate_target(
                    D, N, candidate, target, modular_degree, c_E, old_form,
                    shimura=shimura, homology_matrices=homology_matrices,
                    with_solution=with_solution)
                if result["represents_D"]:
                    return result
            except Exception as error:
                errors.append(_error_record(
                    "target lattice", error, N=N, M=M, E=label,
                    target=target.label, delta=target.degree))
        return None

    positive = first_positive(fixed_targets)
    if positive is not None:
        return [positive], errors
    if not lattice_targets:
        return [], errors

    try:
        shimura = embedded_shimura_data(source_key, N)
    except NotImplementedError as error:
        # Only an unavailable operation triggers the homology fallback.
        warnings.warn("N=%s, E=%s: Shimura shortcut unavailable; "
                      "using integral homology (%s)" % (N, label, error))
        shimura = {"applies": False}
    except Exception as error:
        errors.append(_error_record("Shimura computation", error, N=N, M=M, E=label))
        return [], errors
    homology_matrices = None
    if not shimura["applies"]:
        try:
            homology_divisors, homology_matrices = degeneracy_matrices(source_key, N)
            if tuple(homology_divisors) != tuple(sorted(divisors(N // M))):
                raise ArithmeticError("inconsistent degeneracy-map order")
        except Exception as error:
            errors.append(_error_record("homology setup", error, N=N, M=M, E=label))
            return [], errors
    positive = first_positive(
        lattice_targets, shimura=shimura if shimura["applies"] else None,
        homology_matrices=homology_matrices)
    return ([positive] if positive is not None else []), errors


def _error_record(stage, error, **context):
    """Keep a failed computation distinct from a negative decision."""
    return {
        "stage": stage, "message": str(error),
        **{key: int(value) if isinstance(value, (int, Integer)) else str(value)
           for key, value in context.items()},
    }


# Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("D", type=int, help="positive map degree")
    parser.add_argument("--with-solution", action="store_true", help="also find one integral solution")
    args = parser.parse_args()
    if args.D <= 0:
        parser.error("D must be a positive integer")
    sys.path.insert(0, str(SOURCE_DIR))
    from d_elliptic_output import format_results
    result = compute_degree(args.D, with_solution=args.with_solution)
    print(format_results(result), end="")
    return 0 if result["status"] == "complete" else 1


# Sage CLI and library loads use different entry paths.
if not globals().get("_loaded_as_library", False):
    raise SystemExit(int(main()))
