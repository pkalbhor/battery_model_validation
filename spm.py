import pybamm
import numpy as np
from sei_model import CustomSEI

# pybamm.set_logging_level("DEBUG")

# load lithium-ion base model
model = pybamm.lithium_ion.BaseModel(name="my li-ion model")

# set options
options = {"SEI": "reaction limited", 'lithium plating': 'partially reversible'}
print("Options", model.options)
model.options.update(options)

# set choice of submodels
model.submodels["external circuit"] = pybamm.external_circuit.ExplicitCurrentControl(
    model.param, model.options
)
model.submodels["discharge and throughput variables"] = (
            pybamm.external_circuit.DischargeThroughput(model.param, model.options)
        )
model.submodels["current collector"] = pybamm.current_collector.Uniform(model.param)
model.submodels["thermal"] = pybamm.thermal.isothermal.Isothermal(model.param)
model.submodels["porosity"] = pybamm.porosity.Constant(model.param, model.options)
model.submodels["electrolyte diffusion"] = (
    pybamm.electrolyte_diffusion.ConstantConcentration(model.param)
)
model.submodels["electrolyte conductivity"] = (
    pybamm.electrolyte_conductivity.LeadingOrder(model.param)
)

# Loop over negative and positive electrode domains for some submodels
for domain in ["negative", "positive"]:
    model.submodels[f"{domain} active material"] = pybamm.active_material.Constant(
        model.param, domain, model.options
    )
    model.submodels[f"{domain} electrode potential"] = (
        pybamm.electrode.ohm.LeadingOrder(model.param, domain)
    )
    # model.submodels[f"{domain} particle"] = pybamm.particle.XAveragedPolynomialProfile(
    #     model.param,
    #     domain,
    #     options={**model.options, "particle": "uniform profile"},
    #     phase="primary",
    # )
    model.submodels[f"{domain} particle"] = pybamm.particle.FickianDiffusion(
        model.param,
        domain,
        options={**model.options, "particle": "Fickian diffusion"},
        phase="primary",
    )
    model.submodels[f"{domain} total particle concentration"] = (
        pybamm.particle.TotalConcentration(
            model.param, domain, model.options, phase="primary"
        )
    )

    model.submodels[f"{domain} open-circuit potential"] = (
        pybamm.open_circuit_potential.SingleOpenCircuitPotential(
            model.param,
            domain,
            "lithium-ion main",
            options=model.options,
            phase="primary",
        )
    )
    model.submodels[f"{domain} interface"] = pybamm.kinetics.InverseButlerVolmer(
        model.param, domain, "lithium-ion main", options=model.options
    )
    model.submodels[f"{domain} interface utilisation"] = (
        pybamm.interface_utilisation.Full(model.param, domain, model.options)
    )
    model.submodels[f"{domain} interface current"] = (
        pybamm.kinetics.CurrentForInverseButlerVolmer(
            model.param, domain, "lithium-ion main"
        )
    )
    model.submodels[f"{domain} surface potential difference [V]"] = (
        pybamm.electrolyte_conductivity.surface_potential_form.Explicit(
            model.param, domain, model.options
        )
    )
    model.submodels[f"{domain} particle mechanics"] = (
        pybamm.particle_mechanics.NoMechanics(model.param, domain, model.options)
    )
    model.submodels[f"{domain} sei"] = pybamm.sei.SEIGrowth(
        model.param, domain, "full electrode", model.options, "primary"
    )
    model.submodels[f"{domain} sei on cracks"] = pybamm.sei.NoSEI(
        model.param, domain, model.options, cracks=True
    )
    model.submodels[f"{domain} lithium plating"] = pybamm.lithium_plating.NoPlating(
        model.param, domain, model.options
    )
    model.submodels[f"{domain} lithium plating"] = pybamm.lithium_plating.Plating(
        model.param, domain, True, model.options
    )

# build model
model.build_model()

# solve model
t_eval = np.linspace(0, 3600, 100)

parameter_values = pybamm.ParameterValues("PKalbhor2024")
sim = pybamm.Simulation(model=model, parameter_values=parameter_values)
if False:
    # create geometry
    geometry = pybamm.battery_geometry()

    # process model and geometry
    param = model.default_parameter_values
    param.process_model(model)
    param.process_geometry(geometry)

    # set mesh
    # Note: li-ion base model has defaults for mesh and var_pts
    mesh = pybamm.Mesh(geometry, model.default_submesh_types, model.default_var_pts)

    # discretise model
    # Note: li-ion base model has default spatial methods
    disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
    disc.process_model(model)

    solver = pybamm.ScipySolver()
    solution = solver.solve(model, t_eval)
