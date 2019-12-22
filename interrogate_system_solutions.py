import pybindingcurve as pbc
import random

analytical_system = pbc.BindingCurve(pbc.systems.System_analytical_homodimerbreaking_pp)
kinetic_system = pbc.BindingCurve(pbc.systems.System_kinetic_homodimerbreaking_pp)

query_system={
    'p':random.uniform(0.0,1000.0)*10**random.choice([0,-3,-6,-9,-12]),
    'l':random.uniform(0.0,1000.0)*10**random.choice([0,-3,-6,-9,-12]),
    'kdpl':random.uniform(0.0,1000.0)*10**random.choice([0,-3,-6,-9,-12]),
    'kdpp':random.uniform(0.0,1000.0)*10**random.choice([0,-3,-6,-9,-12])
    }
print(query_system)
kinetic_result=kinetic_system.query(query_system)
analytical_result=analytical_system.query(query_system)

print(f"Kinetic {kinetic_result}")
print(f"Analytical {analytical_result}")

closest = min(range(len(analytical_result)), key=lambda i: abs(analytical_result[i]-kinetic_result))
print(f"Index = {closest}")