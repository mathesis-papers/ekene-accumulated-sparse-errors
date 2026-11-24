import json
from fractions import Fraction


def error_value(dimension, coefficients):
    coefficient_sum = sum(coefficients, Fraction())
    coordinate_error = sum(
        ((Fraction(1, dimension) - value) ** 2 for value in coefficients),
        Fraction(),
    )
    return ((1 - coefficient_sum) ** 2 + coordinate_error) / 2


def optimum(dimension, degree):
    return Fraction(
        (dimension + 1) * (dimension - degree),
        2 * dimension ** 2 * (degree + 1),
    )


def rational(value):
    return f"{value.numerator}/{value.denominator}"


def parameter_error(dimension, correlation, coefficients):
    independent = 1 - correlation
    coefficient_sum = sum(coefficients, Fraction())
    coordinate_error = sum(
        ((Fraction(1, dimension) - value) ** 2 for value in coefficients),
        Fraction(),
    )
    return independent * coordinate_error + correlation * (1 - coefficient_sum) ** 2


def parameter_optimum(dimension, degree, correlation):
    independent = 1 - correlation
    numerator = independent * (independent + dimension * correlation) * (dimension - degree)
    return numerator / (dimension ** 2 * (independent + degree * correlation))


def parameter_checks():
    dimensions = list(range(1, 17)) + [32, 64]
    correlations = [Fraction(1, 17), Fraction(1, 5), Fraction(1, 2), Fraction(4, 5), Fraction(16, 17)]
    optimizer_checks = 0
    identity_checks = 0
    gram_checks = 0
    total_checks = 0
    normal_checks = 0
    for dimension in dimensions:
        for correlation in correlations:
            independent = 1 - correlation
            full_factor = independent + dimension * correlation
            assert 0 < full_factor / dimension <= 1
            shift = independent / correlation
            direct_total = sum((parameter_optimum(dimension, degree, correlation)
                                for degree in range(1, dimension + 1)), Fraction())
            shifted_sum = sum((1 / (degree + shift)
                               for degree in range(1, dimension + 1)), Fraction())
            closed_total = independent * (dimension + shift) / dimension ** 2
            closed_total *= (dimension + shift) * shifted_sum - dimension
            assert direct_total == closed_total
            total_checks += 1
            for degree in range(1, dimension + 1):
                sparse_factor = independent + degree * correlation
                coefficient = full_factor / (dimension * sparse_factor)
                coefficients = [coefficient] * degree + [Fraction()] * (dimension - degree)
                optimum_value = parameter_optimum(dimension, degree, correlation)
                assert parameter_error(dimension, correlation, coefficients) == optimum_value
                assert sparse_factor * coefficient == full_factor / dimension
                assert optimum_value <= Fraction(1, degree)
                optimizer_checks += 1
                for control in [Fraction(-3, 2), Fraction(), Fraction(5, 3)]:
                    candidate = [control * Fraction(index + 1, degree) for index in range(degree)]
                    candidate += [Fraction()] * (dimension - degree)
                    candidate_sum = sum(candidate, Fraction())
                    squares = sum((value ** 2 for value in candidate), Fraction())
                    center = degree * coefficient
                    residual = independent * (squares - candidate_sum ** 2 / degree)
                    residual += (independent / degree + correlation) * (candidate_sum - center) ** 2
                    direct_error = parameter_error(dimension, correlation, candidate)
                    assert direct_error - optimum_value == residual >= 0
                    identity_checks += 1
                    if dimension <= 8:
                        difference = [Fraction(1, dimension) - value for value in candidate]
                        gram_error = sum((difference[row] * difference[column]
                                          * (correlation + independent * (row == column))
                                          for row in range(dimension) for column in range(dimension)), Fraction())
                        assert direct_error == gram_error
                        gram_checks += 1
                if dimension <= 8:
                    augmented = [
                        [correlation + independent * (row == column) for column in range(degree)]
                        + [full_factor / dimension] for row in range(degree)
                    ]
                    for pivot in range(degree):
                        divisor = augmented[pivot][pivot]
                        assert divisor > 0
                        augmented[pivot] = [value / divisor for value in augmented[pivot]]
                        for row in range(degree):
                            if row != pivot:
                                multiplier = augmented[row][pivot]
                                augmented[row] = [value - multiplier * pivot_value
                                                  for value, pivot_value in zip(augmented[row], augmented[pivot])]
                    assert [row[-1] for row in augmented] == [coefficient] * degree
                    normal_checks += 1
                if correlation == Fraction(1, 2):
                    assert optimum_value == optimum(dimension, degree)
            assert parameter_error(dimension, correlation, [Fraction(1, dimension)] * dimension) == 0
    harmonic_checks = []
    for denominator in [2, 3, 4, 8, 16]:
        dimension = 2 ** denominator
        correlation = Fraction(1, denominator)
        shifted_harmonic = sum((Fraction(1, degree + denominator - 1)
                                for degree in range(1, dimension + 1)), Fraction())
        closed_total = Fraction((denominator - 1) * (dimension + denominator - 1), denominator * dimension ** 2)
        closed_total *= (dimension + denominator - 1) * shifted_harmonic - dimension
        if dimension <= 256:
            assert closed_total == sum((parameter_optimum(dimension, degree, correlation)
                                        for degree in range(1, dimension + 1)), Fraction())
        harmonic_checks.append({'dimension': dimension, 'correlation': rational(correlation),
                                'total_lower_floor': closed_total.numerator // closed_total.denominator})
    return {
        'dimensions': dimensions,
        'correlations': [rational(value) for value in correlations],
        'optimizer_and_stationarity_checks': optimizer_checks,
        'signed_square_completion_checks': identity_checks,
        'direct_gram_identity_checks': gram_checks,
        'independent_normal_equation_checks': normal_checks,
        'closed_total_checks': total_checks,
        'shifted_harmonic_cases': harmonic_checks,
        'scope': 'Finite rational parameter checks; the complete parameter range and sharp asymptotic coefficient have written proofs in note.tex.',
    }


def scalar_sampling_checks():
    checks = 0
    for root_degree in [1, 2, 3, 4, 8]:
        degree = root_degree ** 2
        for target_energy, variance in [
            (Fraction(), Fraction()),
            (Fraction(), Fraction(1)),
            (Fraction(1), Fraction()),
            (Fraction(1, 4), Fraction(3, 4)),
            (Fraction(1, 3), Fraction(1, 6)),
            (Fraction(1, root_degree + 1), Fraction(root_degree, root_degree + 1)),
        ]:
            denominator = variance + degree * target_energy
            if denominator:
                scale = degree * target_energy / denominator
                expected_error = (1 - scale) ** 2 * target_energy + scale ** 2 * variance / degree
                assert expected_error == target_energy * variance / denominator
            else:
                expected_error = Fraction()
            assert expected_error <= (target_energy + variance) / (1 + root_degree) ** 2
            checks += 1
    return {'squared_degrees': [1, 4, 9, 16, 64], 'exact_scalar_minimization_cases': checks,
            'scope': 'Rational controls of the scalar expectation and its envelope, including zero-energy and zero-variance endpoints.'}


def solve_normal_equations(dimension, degree):
    augmented = [
        [Fraction(1 + (row == column)) for column in range(degree)]
        + [1 + Fraction(1, dimension)]
        for row in range(degree)
    ]
    for pivot in range(degree):
        divisor = augmented[pivot][pivot]
        assert divisor != 0
        augmented[pivot] = [value / divisor for value in augmented[pivot]]
        for row in range(degree):
            if row == pivot:
                continue
            multiplier = augmented[row][pivot]
            if multiplier:
                augmented[row] = [
                    value - multiplier * pivot_value
                    for value, pivot_value in zip(augmented[row], augmented[pivot])
                ]
    return [row[-1] for row in augmented]


def check():
    dimensions = list(range(1, 25)) + [32, 64, 128]
    optimizer_checks = 0
    identity_checks = 0
    signed_controls = [Fraction(-3, 2), Fraction(-1, 7), Fraction(), Fraction(2, 3), Fraction(5, 2)]
    for dimension in dimensions:
        harmonic = sum((Fraction(1, index) for index in range(1, dimension + 1)), Fraction())
        direct_total = sum((optimum(dimension, degree) for degree in range(1, dimension)), Fraction())
        closed_total = Fraction(dimension + 1, 2 * dimension ** 2) * ((dimension + 1) * harmonic - 2 * dimension)
        assert direct_total == closed_total
        assert Fraction(dimension + 1, 2 * dimension) <= 1
        for degree in range(1, dimension + 1):
            coefficient = Fraction(dimension + 1, dimension * (degree + 1))
            coefficients = [coefficient] * degree + [Fraction()] * (dimension - degree)
            value = error_value(dimension, coefficients)
            assert value == optimum(dimension, degree)
            assert value <= Fraction(1, degree)
            coefficient_sum = degree * coefficient
            assert coefficient + coefficient_sum == 1 + Fraction(1, dimension)
            optimizer_checks += 1
            for control in signed_controls:
                candidate = [control * Fraction(index + 1, degree) for index in range(degree)]
                candidate += [Fraction()] * (dimension - degree)
                candidate_sum = sum(candidate, Fraction())
                squares = sum((value ** 2 for value in candidate), Fraction())
                center = Fraction(degree * (dimension + 1), dimension * (degree + 1))
                residual = (squares - candidate_sum ** 2 / degree) / 2
                residual += Fraction(degree + 1, 2 * degree) * (candidate_sum - center) ** 2
                assert error_value(dimension, candidate) - value == residual
                assert residual >= 0
                identity_checks += 1
        assert error_value(dimension, [Fraction(1, dimension)] * dimension) == 0
    total = sum((optimum(64, degree) for degree in range(1, 64)), Fraction())
    harmonic_64 = sum((Fraction(1, index) for index in range(1, 65)), Fraction())
    assert harmonic_64 >= 4
    assert total >= Fraction(2145, 2048) > 1
    assert Fraction(1431, 1000) < total < Fraction(1432, 1000)
    normal_equation_cases = [(1, 1), (2, 1), (8, 3), (16, 7)]
    normal_equation_cases += [(64, degree) for degree in [1, 2, 8, 32, 64]]
    for dimension, degree in normal_equation_cases:
        solution = solve_normal_equations(dimension, degree)
        assert solution == [Fraction(dimension + 1, dimension * (degree + 1))] * degree
        solution += [Fraction()] * (dimension - degree)
        assert error_value(dimension, solution) == optimum(dimension, degree)
    return {
        "arithmetic": "exact rational arithmetic; no floating-point premises",
        "dimensions": dimensions,
        "optimizer_and_stationarity_checks": optimizer_checks,
        "signed_square_completion_checks": identity_checks,
        "full_atom_representation_checks": len(dimensions),
        "independent_normal_equation_cases": normal_equation_cases,
        "dimension_64_total": rational(total),
        "dimension_64_interval": {"strict_lower": "1431/1000", "strict_upper": "1432/1000"},
        "dyadic_lower_bound": "2145/2048",
        "parameter_family": parameter_checks(),
        "scalar_sampling": scalar_sampling_checks(),
        "coverage": "Finite rational identities and optimizer values; the all-parameter minimum and support equivalence have written proofs in note.tex.",
        "status": "all checks passed",
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
