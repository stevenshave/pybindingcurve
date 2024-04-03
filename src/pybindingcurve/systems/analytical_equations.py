from mpmath import mpf, sqrt, power, fabs, almosteq

# 1:1 binding - see https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
# Readout is PL
def system01_analytical_one_to_one__pl(p, l, kdpl):
    p = mpf(p)
    l = mpf(l)
    kdpl = mpf(kdpl)
    return ((p + kdpl + l - sqrt(-4 * p * l + power(p + kdpl + l, 2))) / 2.0).real


# 1:1:1 competition - see https://stevenshave.github.io/pybindingcurve/simulate_competition.html
# Readout is PL
def system02_analytical_competition__pl(p, l, i, kdpl, kdpi):
    p = mpf(p)
    l = mpf(l)
    i = mpf(i)
    kdpl = mpf(kdpl)
    kdpi = mpf(kdpi)
    if almosteq(kdpl, kdpi, 1e-10):
        kdpi += 1e-10
    if kdpl < kdpi:
        return (
            -(
                p * kdpi
                - p * kdpl
                + i * kdpl
                + kdpi * kdpl
                - power(kdpl, 2)
                + 2 * kdpi * l
                - kdpl * l
            )
            / (3.0 * (-kdpi + kdpl))
            - (
                power(2, 0.3333333333333333)
                * (
                    -power(
                        p * kdpi
                        - p * kdpl
                        + i * kdpl
                        + kdpi * kdpl
                        - power(kdpl, 2)
                        + 2 * kdpi * l
                        - kdpl * l,
                        2,
                    )
                    + 3
                    * (-kdpi + kdpl)
                    * (
                        -2 * p * kdpi * l
                        + p * kdpl * l
                        - i * kdpl * l
                        - kdpi * kdpl * l
                        - kdpi * power(l, 2)
                    )
                )
            )
            / (
                3.0
                * (-kdpi + kdpl)
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            + power(
                -2 * power(p, 3) * power(kdpi, 3)
                + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                + 2 * power(p, 3) * power(kdpl, 3)
                - 6 * power(p, 2) * i * power(kdpl, 3)
                + 6 * p * power(i, 2) * power(kdpl, 3)
                - 2 * power(i, 3) * power(kdpl, 3)
                - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                + 24 * p * i * kdpi * power(kdpl, 3)
                - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                - 2 * power(kdpi, 3) * power(kdpl, 3)
                + 6 * power(p, 2) * power(kdpl, 4)
                - 12 * p * i * power(kdpl, 4)
                + 6 * power(i, 2) * power(kdpl, 4)
                - 18 * p * kdpi * power(kdpl, 4)
                + 12 * i * kdpi * power(kdpl, 4)
                + 6 * power(kdpi, 2) * power(kdpl, 4)
                + 6 * p * power(kdpl, 5)
                - 6 * i * power(kdpl, 5)
                - 6 * kdpi * power(kdpl, 5)
                + 2 * power(kdpl, 6)
                + 6 * power(p, 2) * power(kdpi, 3) * l
                - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                + 3 * p * i * power(kdpi, 2) * kdpl * l
                + 3 * p * power(kdpi, 3) * kdpl * l
                + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                - 9 * p * i * kdpi * power(kdpl, 2) * l
                - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                - 3 * power(p, 2) * power(kdpl, 3) * l
                + 6 * p * i * power(kdpl, 3) * l
                - 3 * power(i, 2) * power(kdpl, 3) * l
                - 3 * p * kdpi * power(kdpl, 3) * l
                + 9 * i * kdpi * power(kdpl, 3) * l
                + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                + 3 * p * power(kdpl, 4) * l
                - 3 * i * power(kdpl, 4) * l
                - 15 * kdpi * power(kdpl, 4) * l
                + 6 * power(kdpl, 5) * l
                - 6 * p * power(kdpi, 3) * power(l, 2)
                + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                - 3 * p * power(kdpl, 3) * power(l, 2)
                + 3 * i * power(kdpl, 3) * power(l, 2)
                - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                + 6 * power(kdpl, 4) * power(l, 2)
                + 2 * power(kdpi, 3) * power(l, 3)
                - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                + 2 * power(kdpl, 3) * power(l, 3)
                + sqrt(
                    power(
                        -2 * power(p, 3) * power(kdpi, 3)
                        + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                        - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                        - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                        - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                        + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                        - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                        + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                        - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                        - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                        + 2 * power(p, 3) * power(kdpl, 3)
                        - 6 * power(p, 2) * i * power(kdpl, 3)
                        + 6 * p * power(i, 2) * power(kdpl, 3)
                        - 2 * power(i, 3) * power(kdpl, 3)
                        - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                        + 24 * p * i * kdpi * power(kdpl, 3)
                        - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                        + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                        - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                        - 2 * power(kdpi, 3) * power(kdpl, 3)
                        + 6 * power(p, 2) * power(kdpl, 4)
                        - 12 * p * i * power(kdpl, 4)
                        + 6 * power(i, 2) * power(kdpl, 4)
                        - 18 * p * kdpi * power(kdpl, 4)
                        + 12 * i * kdpi * power(kdpl, 4)
                        + 6 * power(kdpi, 2) * power(kdpl, 4)
                        + 6 * p * power(kdpl, 5)
                        - 6 * i * power(kdpl, 5)
                        - 6 * kdpi * power(kdpl, 5)
                        + 2 * power(kdpl, 6)
                        + 6 * power(p, 2) * power(kdpi, 3) * l
                        - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                        + 3 * p * i * power(kdpi, 2) * kdpl * l
                        + 3 * p * power(kdpi, 3) * kdpl * l
                        + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                        - 9 * p * i * kdpi * power(kdpl, 2) * l
                        - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                        - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                        - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                        - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                        - 3 * power(p, 2) * power(kdpl, 3) * l
                        + 6 * p * i * power(kdpl, 3) * l
                        - 3 * power(i, 2) * power(kdpl, 3) * l
                        - 3 * p * kdpi * power(kdpl, 3) * l
                        + 9 * i * kdpi * power(kdpl, 3) * l
                        + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                        + 3 * p * power(kdpl, 4) * l
                        - 3 * i * power(kdpl, 4) * l
                        - 15 * kdpi * power(kdpl, 4) * l
                        + 6 * power(kdpl, 5) * l
                        - 6 * p * power(kdpi, 3) * power(l, 2)
                        + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                        + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                        + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                        - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                        - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                        + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                        - 3 * p * power(kdpl, 3) * power(l, 2)
                        + 3 * i * power(kdpl, 3) * power(l, 2)
                        - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                        + 6 * power(kdpl, 4) * power(l, 2)
                        + 2 * power(kdpi, 3) * power(l, 3)
                        - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                        - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                        + 2 * power(kdpl, 3) * power(l, 3),
                        2,
                    )
                    + 4
                    * power(
                        -power(
                            p * kdpi
                            - p * kdpl
                            + i * kdpl
                            + kdpi * kdpl
                            - power(kdpl, 2)
                            + 2 * kdpi * l
                            - kdpl * l,
                            2,
                        )
                        + 3
                        * (-kdpi + kdpl)
                        * (
                            -2 * p * kdpi * l
                            + p * kdpl * l
                            - i * kdpl * l
                            - kdpi * kdpl * l
                            - kdpi * power(l, 2)
                        ),
                        3,
                    )
                ),
                0.3333333333333333,
            )
            / (3.0 * power(2, 0.3333333333333333) * (-kdpi + kdpl))
        ).real
    else:
        return (
            -(
                p * kdpi
                - p * kdpl
                + i * kdpl
                + kdpi * kdpl
                - power(kdpl, 2)
                + 2 * kdpi * l
                - kdpl * l
            )
            / (3.0 * (-kdpi + kdpl))
            + (
                (1 - complex(0, 1) * sqrt(3))
                * (
                    -power(
                        p * kdpi
                        - p * kdpl
                        + i * kdpl
                        + kdpi * kdpl
                        - power(kdpl, 2)
                        + 2 * kdpi * l
                        - kdpl * l,
                        2,
                    )
                    + 3
                    * (-kdpi + kdpl)
                    * (
                        -2 * p * kdpi * l
                        + p * kdpl * l
                        - i * kdpl * l
                        - kdpi * kdpl * l
                        - kdpi * power(l, 2)
                    )
                )
            )
            / (
                3.0
                * power(2, 0.6666666666666666)
                * (-kdpi + kdpl)
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            - (
                (1 + complex(0, 1) * sqrt(3))
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            / (6.0 * power(2, 0.3333333333333333) * (-kdpi + kdpl))
        ).real


# Homodimer formation - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
# Readout is PP
def system03_analytical_homodimer_formation__pp(p, kdpp):
    p = mpf(p)
    kdpp = mpf(kdpp)
    return ((4 * p + kdpp - sqrt(kdpp) * sqrt(8 * p + kdpp)) / 8.0).real


# Homodimer breaking - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
# Readout is PP
def system04_analytical_homodimer_breaking__pp(p, i, kdpp, kdpi):
    p = mpf(p)
    i = mpf(i)
    kdpp = mpf(kdpp)
    kdpi = mpf(kdpi)
    return [
        (
            -(-4 * p * kdpp - power(kdpp, 2) - 4 * power(kdpi, 2) + 4 * kdpp * i)
            / (12.0 * kdpp)
            - (
                -power(
                    -4 * p * kdpp - power(kdpp, 2) - 4 * power(kdpi, 2) + 4 * kdpp * i,
                    2,
                )
                + 12
                * kdpp
                * (
                    power(p, 2) * kdpp
                    + 4 * p * power(kdpi, 2)
                    + kdpp * power(kdpi, 2)
                    - 2 * p * kdpp * i
                    + 2 * kdpp * kdpi * i
                    + kdpp * power(i, 2)
                )
            )
            / (
                6.0
                * power(2, 0.6666666666666666)
                * kdpp
                * power(
                    -16 * power(p, 3) * power(kdpp, 3)
                    + 60 * power(p, 2) * power(kdpp, 4)
                    + 24 * p * power(kdpp, 5)
                    + 2 * power(kdpp, 6)
                    + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                    - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                    - 12 * power(kdpp, 4) * power(kdpi, 2)
                    - 192 * p * kdpp * power(kdpi, 4)
                    - 48 * power(kdpp, 2) * power(kdpi, 4)
                    + 128 * power(kdpi, 6)
                    + 48 * power(p, 2) * power(kdpp, 3) * i
                    - 120 * p * power(kdpp, 4) * i
                    - 24 * power(kdpp, 5) * i
                    - 288 * p * power(kdpp, 3) * kdpi * i
                    - 72 * power(kdpp, 4) * kdpi * i
                    + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                    - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                    - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                    - 384 * kdpp * power(kdpi, 4) * i
                    - 48 * p * power(kdpp, 3) * power(i, 2)
                    + 60 * power(kdpp, 4) * power(i, 2)
                    + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                    + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                    + 16 * power(kdpp, 3) * power(i, 3)
                    + sqrt(
                        power(
                            -16 * power(p, 3) * power(kdpp, 3)
                            + 60 * power(p, 2) * power(kdpp, 4)
                            + 24 * p * power(kdpp, 5)
                            + 2 * power(kdpp, 6)
                            + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                            - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                            - 12 * power(kdpp, 4) * power(kdpi, 2)
                            - 192 * p * kdpp * power(kdpi, 4)
                            - 48 * power(kdpp, 2) * power(kdpi, 4)
                            + 128 * power(kdpi, 6)
                            + 48 * power(p, 2) * power(kdpp, 3) * i
                            - 120 * p * power(kdpp, 4) * i
                            - 24 * power(kdpp, 5) * i
                            - 288 * p * power(kdpp, 3) * kdpi * i
                            - 72 * power(kdpp, 4) * kdpi * i
                            + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                            - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                            - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                            - 384 * kdpp * power(kdpi, 4) * i
                            - 48 * p * power(kdpp, 3) * power(i, 2)
                            + 60 * power(kdpp, 4) * power(i, 2)
                            + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                            + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                            + 16 * power(kdpp, 3) * power(i, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                -4 * p * kdpp
                                - power(kdpp, 2)
                                - 4 * power(kdpi, 2)
                                + 4 * kdpp * i,
                                2,
                            )
                            + 12
                            * kdpp
                            * (
                                power(p, 2) * kdpp
                                + 4 * p * power(kdpi, 2)
                                + kdpp * power(kdpi, 2)
                                - 2 * p * kdpp * i
                                + 2 * kdpp * kdpi * i
                                + kdpp * power(i, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            + power(
                -16 * power(p, 3) * power(kdpp, 3)
                + 60 * power(p, 2) * power(kdpp, 4)
                + 24 * p * power(kdpp, 5)
                + 2 * power(kdpp, 6)
                + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                - 12 * power(kdpp, 4) * power(kdpi, 2)
                - 192 * p * kdpp * power(kdpi, 4)
                - 48 * power(kdpp, 2) * power(kdpi, 4)
                + 128 * power(kdpi, 6)
                + 48 * power(p, 2) * power(kdpp, 3) * i
                - 120 * p * power(kdpp, 4) * i
                - 24 * power(kdpp, 5) * i
                - 288 * p * power(kdpp, 3) * kdpi * i
                - 72 * power(kdpp, 4) * kdpi * i
                + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                - 384 * kdpp * power(kdpi, 4) * i
                - 48 * p * power(kdpp, 3) * power(i, 2)
                + 60 * power(kdpp, 4) * power(i, 2)
                + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                + 16 * power(kdpp, 3) * power(i, 3)
                + sqrt(
                    power(
                        -16 * power(p, 3) * power(kdpp, 3)
                        + 60 * power(p, 2) * power(kdpp, 4)
                        + 24 * p * power(kdpp, 5)
                        + 2 * power(kdpp, 6)
                        + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                        - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                        - 12 * power(kdpp, 4) * power(kdpi, 2)
                        - 192 * p * kdpp * power(kdpi, 4)
                        - 48 * power(kdpp, 2) * power(kdpi, 4)
                        + 128 * power(kdpi, 6)
                        + 48 * power(p, 2) * power(kdpp, 3) * i
                        - 120 * p * power(kdpp, 4) * i
                        - 24 * power(kdpp, 5) * i
                        - 288 * p * power(kdpp, 3) * kdpi * i
                        - 72 * power(kdpp, 4) * kdpi * i
                        + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                        - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                        - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                        - 384 * kdpp * power(kdpi, 4) * i
                        - 48 * p * power(kdpp, 3) * power(i, 2)
                        + 60 * power(kdpp, 4) * power(i, 2)
                        + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                        + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                        + 16 * power(kdpp, 3) * power(i, 3),
                        2,
                    )
                    + 4
                    * power(
                        -power(
                            -4 * p * kdpp
                            - power(kdpp, 2)
                            - 4 * power(kdpi, 2)
                            + 4 * kdpp * i,
                            2,
                        )
                        + 12
                        * kdpp
                        * (
                            power(p, 2) * kdpp
                            + 4 * p * power(kdpi, 2)
                            + kdpp * power(kdpi, 2)
                            - 2 * p * kdpp * i
                            + 2 * kdpp * kdpi * i
                            + kdpp * power(i, 2)
                        ),
                        3,
                    )
                ),
                0.3333333333333333,
            )
            / (12.0 * power(2, 0.3333333333333333) * kdpp)
        ).real,
        (
            -(-4 * p * kdpp - power(kdpp, 2) - 4 * power(kdpi, 2) + 4 * kdpp * i)
            / (12.0 * kdpp)
            + (
                (1 - complex(0, 1) * sqrt(3))
                * (
                    -power(
                        -4 * p * kdpp
                        - power(kdpp, 2)
                        - 4 * power(kdpi, 2)
                        + 4 * kdpp * i,
                        2,
                    )
                    + 12
                    * kdpp
                    * (
                        power(p, 2) * kdpp
                        + 4 * p * power(kdpi, 2)
                        + kdpp * power(kdpi, 2)
                        - 2 * p * kdpp * i
                        + 2 * kdpp * kdpi * i
                        + kdpp * power(i, 2)
                    )
                )
            )
            / (
                12.0
                * power(2, 0.6666666666666666)
                * kdpp
                * power(
                    -16 * power(p, 3) * power(kdpp, 3)
                    + 60 * power(p, 2) * power(kdpp, 4)
                    + 24 * p * power(kdpp, 5)
                    + 2 * power(kdpp, 6)
                    + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                    - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                    - 12 * power(kdpp, 4) * power(kdpi, 2)
                    - 192 * p * kdpp * power(kdpi, 4)
                    - 48 * power(kdpp, 2) * power(kdpi, 4)
                    + 128 * power(kdpi, 6)
                    + 48 * power(p, 2) * power(kdpp, 3) * i
                    - 120 * p * power(kdpp, 4) * i
                    - 24 * power(kdpp, 5) * i
                    - 288 * p * power(kdpp, 3) * kdpi * i
                    - 72 * power(kdpp, 4) * kdpi * i
                    + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                    - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                    - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                    - 384 * kdpp * power(kdpi, 4) * i
                    - 48 * p * power(kdpp, 3) * power(i, 2)
                    + 60 * power(kdpp, 4) * power(i, 2)
                    + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                    + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                    + 16 * power(kdpp, 3) * power(i, 3)
                    + sqrt(
                        power(
                            -16 * power(p, 3) * power(kdpp, 3)
                            + 60 * power(p, 2) * power(kdpp, 4)
                            + 24 * p * power(kdpp, 5)
                            + 2 * power(kdpp, 6)
                            + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                            - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                            - 12 * power(kdpp, 4) * power(kdpi, 2)
                            - 192 * p * kdpp * power(kdpi, 4)
                            - 48 * power(kdpp, 2) * power(kdpi, 4)
                            + 128 * power(kdpi, 6)
                            + 48 * power(p, 2) * power(kdpp, 3) * i
                            - 120 * p * power(kdpp, 4) * i
                            - 24 * power(kdpp, 5) * i
                            - 288 * p * power(kdpp, 3) * kdpi * i
                            - 72 * power(kdpp, 4) * kdpi * i
                            + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                            - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                            - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                            - 384 * kdpp * power(kdpi, 4) * i
                            - 48 * p * power(kdpp, 3) * power(i, 2)
                            + 60 * power(kdpp, 4) * power(i, 2)
                            + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                            + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                            + 16 * power(kdpp, 3) * power(i, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                -4 * p * kdpp
                                - power(kdpp, 2)
                                - 4 * power(kdpi, 2)
                                + 4 * kdpp * i,
                                2,
                            )
                            + 12
                            * kdpp
                            * (
                                power(p, 2) * kdpp
                                + 4 * p * power(kdpi, 2)
                                + kdpp * power(kdpi, 2)
                                - 2 * p * kdpp * i
                                + 2 * kdpp * kdpi * i
                                + kdpp * power(i, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            - (
                (1 + complex(0, 1) * sqrt(3))
                * power(
                    -16 * power(p, 3) * power(kdpp, 3)
                    + 60 * power(p, 2) * power(kdpp, 4)
                    + 24 * p * power(kdpp, 5)
                    + 2 * power(kdpp, 6)
                    + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                    - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                    - 12 * power(kdpp, 4) * power(kdpi, 2)
                    - 192 * p * kdpp * power(kdpi, 4)
                    - 48 * power(kdpp, 2) * power(kdpi, 4)
                    + 128 * power(kdpi, 6)
                    + 48 * power(p, 2) * power(kdpp, 3) * i
                    - 120 * p * power(kdpp, 4) * i
                    - 24 * power(kdpp, 5) * i
                    - 288 * p * power(kdpp, 3) * kdpi * i
                    - 72 * power(kdpp, 4) * kdpi * i
                    + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                    - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                    - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                    - 384 * kdpp * power(kdpi, 4) * i
                    - 48 * p * power(kdpp, 3) * power(i, 2)
                    + 60 * power(kdpp, 4) * power(i, 2)
                    + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                    + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                    + 16 * power(kdpp, 3) * power(i, 3)
                    + sqrt(
                        power(
                            -16 * power(p, 3) * power(kdpp, 3)
                            + 60 * power(p, 2) * power(kdpp, 4)
                            + 24 * p * power(kdpp, 5)
                            + 2 * power(kdpp, 6)
                            + 96 * power(p, 2) * power(kdpp, 2) * power(kdpi, 2)
                            - 96 * p * power(kdpp, 3) * power(kdpi, 2)
                            - 12 * power(kdpp, 4) * power(kdpi, 2)
                            - 192 * p * kdpp * power(kdpi, 4)
                            - 48 * power(kdpp, 2) * power(kdpi, 4)
                            + 128 * power(kdpi, 6)
                            + 48 * power(p, 2) * power(kdpp, 3) * i
                            - 120 * p * power(kdpp, 4) * i
                            - 24 * power(kdpp, 5) * i
                            - 288 * p * power(kdpp, 3) * kdpi * i
                            - 72 * power(kdpp, 4) * kdpi * i
                            + 96 * p * power(kdpp, 2) * power(kdpi, 2) * i
                            - 48 * power(kdpp, 3) * power(kdpi, 2) * i
                            - 288 * power(kdpp, 2) * power(kdpi, 3) * i
                            - 384 * kdpp * power(kdpi, 4) * i
                            - 48 * p * power(kdpp, 3) * power(i, 2)
                            + 60 * power(kdpp, 4) * power(i, 2)
                            + 288 * power(kdpp, 3) * kdpi * power(i, 2)
                            + 240 * power(kdpp, 2) * power(kdpi, 2) * power(i, 2)
                            + 16 * power(kdpp, 3) * power(i, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                -4 * p * kdpp
                                - power(kdpp, 2)
                                - 4 * power(kdpi, 2)
                                + 4 * kdpp * i,
                                2,
                            )
                            + 12
                            * kdpp
                            * (
                                power(p, 2) * kdpp
                                + 4 * p * power(kdpi, 2)
                                + kdpp * power(kdpi, 2)
                                - 2 * p * kdpp * i
                                + 2 * kdpp * kdpi * i
                                + kdpp * power(i, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            / (24.0 * power(2, 0.3333333333333333) * kdpp)
        ).real,
    ]


## Homodimerbreaking-analytical
# def all_solutions(a,x,kdaa,kdax):
# 	return np.array([
# 	-(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x)/(12.*kdaa)-(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)))/(6.*np.power(2,0.6666666666666666)*kdaa*np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333))+np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333)/(12.*np.power(2,0.3333333333333333)*kdaa),
# 	-(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x)/(12.*kdaa)+((1+complex(0,1)*np.lib.scimath.sqrt(3))*(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2))))/(12.*np.power(2,0.6666666666666666)*kdaa*np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333))-((1-complex(0,1)*np.lib.scimath.sqrt(3))*np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333))/(24.*np.power(2,0.3333333333333333)*kdaa),
# 	-(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x)/(12.*kdaa)+((1-complex(0,1)*np.lib.scimath.sqrt(3))*(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2))))/(12.*np.power(2,0.6666666666666666)*kdaa*np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333))-((1+complex(0,1)*np.lib.scimath.sqrt(3))*np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3)+np.lib.scimath.sqrt(np.power(-16*np.power(a,3)*np.power(kdaa,3)+60*np.power(a,2)*np.power(kdaa,4)+24*a*np.power(kdaa,5)+2*np.power(kdaa,6)+96*np.power(a,2)*np.power(kdaa,2)*np.power(kdax,2)-96*a*np.power(kdaa,3)*np.power(kdax,2)-12*np.power(kdaa,4)*np.power(kdax,2)-192*a*kdaa*np.power(kdax,4)-48*np.power(kdaa,2)*np.power(kdax,4)+128*np.power(kdax,6)+48*np.power(a,2)*np.power(kdaa,3)*x-120*a*np.power(kdaa,4)*x-24*np.power(kdaa,5)*x-288*a*np.power(kdaa,3)*kdax*x-72*np.power(kdaa,4)*kdax*x+96*a*np.power(kdaa,2)*np.power(kdax,2)*x-48*np.power(kdaa,3)*np.power(kdax,2)*x-288*np.power(kdaa,2)*np.power(kdax,3)*x-384*kdaa*np.power(kdax,4)*x-48*a*np.power(kdaa,3)*np.power(x,2)+60*np.power(kdaa,4)*np.power(x,2)+288*np.power(kdaa,3)*kdax*np.power(x,2)+240*np.power(kdaa,2)*np.power(kdax,2)*np.power(x,2)+16*np.power(kdaa,3)*np.power(x,3),2)+4*np.power(-np.power(-4*a*kdaa-np.power(kdaa,2)-4*np.power(kdax,2)+4*kdaa*x,2)+12*kdaa*(np.power(a,2)*kdaa+4*a*np.power(kdax,2)+kdaa*np.power(kdax,2)-2*a*kdaa*x+2*kdaa*kdax*x+kdaa*np.power(x,2)),3)),0.3333333333333333))/(24.*np.power(2,0.3333333333333333)*kdaa),
# 	])
