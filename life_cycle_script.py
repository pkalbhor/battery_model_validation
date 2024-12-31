import pandas as pd
import pybamm
from spm import sim, model

save_at_cycles = 20
number_of_cycles = 2000
experiment = pybamm.Experiment(
    operating_conditions=[
        ("Discharge at 1C until 2.0 V",
        "Rest for 2 hour",
        "Discharge at 0.5C until 2.0 V",
        "Rest for 1 hour",
        "Charge at 0.5C until 3.6 V",
        "Rest for 1 hour"),
    ]*number_of_cycles,
    period="60 s"
)

parameters = pybamm.ParameterValues('PKalbhor2024')

sim = pybamm.Simulation(model=model, parameter_values=parameters, experiment=experiment)
sol = sim.solve(initial_soc="3.6 V", showprogress=True, save_at_cycles=save_at_cycles)

cycle_index = list()
max_capacity = list()
for index, cycle in enumerate(sol.cycles):
    if not cycle: continue
    capacity = max(cycle['Discharge capacity [A.h]'].entries)
    max_capacity.append(capacity)
    cycle_index.append(index+1)

df = pd.DataFrame({
    "cycle": cycle_index,
    "max_capacity": max_capacity
})

df.to_csv("live_cycle_data.csv", index=False)
