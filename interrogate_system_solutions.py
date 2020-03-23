import pybindingcurve as pbc
import random

output_file = open("system_results.csv", "w", buffering=1)
output_file.write("p,l,kdpp,kdpl,solution")
analytical_system = pbc.BindingCurve(pbc.systems.System_analytical_homodimerbreaking_pp)
kinetic_system = pbc.BindingCurve(pbc.systems.System_kinetic_homodimerbreaking_pp)

while True:
    query_system = {
        "p": random.uniform(0.0, 1000.0) * 10 ** random.choice([0, -3, -6, -9, -12]),
        "l": random.uniform(0.0, 1000.0) * 10 ** random.choice([0, -3, -6, -9, -12]),
        "kdpl": random.uniform(0.0, 1000.0) * 10 ** random.choice([0, -3, -6, -9, -12]),
        "kdpp": random.uniform(0.0, 1000.0) * 10 ** random.choice([0, -3, -6, -9, -12]),
    }
    kinetic_result = kinetic_system.query(query_system)
    analytical_result = analytical_system.query(query_system)

    closest = min(
        range(len(analytical_result)),
        key=lambda i: abs(analytical_result[i] - kinetic_result),
    )
    output_file.write(
        f"{query_system['p']},{query_system['l']},{query_system['kdpp']},{query_system['kdpl']},{closest}\n"
    )

output_file.close()
