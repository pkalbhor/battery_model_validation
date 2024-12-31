import pybamm

def D(cc):
    return pybamm.FunctionParameter(
        "Diffusivity [m2.s-1]", {"Solvent concentration [mol.m-3]": cc}
    )

class CustomSEI(pybamm.BaseSubModel):
    """
    """
    def __init__(self, param, domain, options=None):
        super().__init__(param, domain, options=options)
        
    def get_fundamental_variables(self):
        xi = pybamm.SpatialVariable("xi", domain="SEI layer", coord_sys="cartesian")
        c = pybamm.Variable("Solvent concentration [mol.m-3]", domain="SEI layer")
        L = pybamm.Variable("SEI thickness [m]")
        
        variables = {
            "xi": xi,
            "Solvent concentration [mol.m-3]": c,
            "SEI thickness [m]": L
        }
        return variables
        
    def get_coupled_variables(self, variables):
        return variables

    def set_rhs(self, variables):

        # dimensional parameters
        k = pybamm.Parameter("Reaction rate constant [m.s-1]")
        V_hat = pybamm.Parameter("Partial molar volume [m3.mol-1]")

        c = variables["Solvent concentration [mol.m-3]"]
        L = variables["SEI thickness [m]"]
        xi = variables["xi"]

        # SEI reaction flux
        R = k * pybamm.BoundaryValue(c, "left")

        # solvent concentration equation
        N = -1 / L * D(c) * pybamm.grad(c)
        dcdt = (V_hat * R) / L * pybamm.inner(xi, pybamm.grad(c)) - 1 / L * pybamm.div(N)

        # SEI thickness equation
        dLdt = V_hat * R
    
        self.rhs = {c: dcdt, L: dLdt}
    

    def set_boundary_conditions(self, variables):
        k = pybamm.Parameter("Reaction rate constant [m.s-1]")
        c_inf = pybamm.Parameter("Bulk electrolyte solvent concentration [mol.m-3]")

        c = variables["Solvent concentration [mol.m-3]"]
        L = variables["SEI thickness [m]"]
        D_left = pybamm.BoundaryValue(
            D(c), "left"
        )  # pybamm requires BoundaryValue(D(c)) and not D(BoundaryValue(c))
        grad_c_left = R * L / D_left
        R = k * pybamm.BoundaryValue(c, "left")
        
        c_right = c_inf
        self.boundary_conditions = {
            c: {"left": (grad_c_left, "Neumann"), "right": (c_right, "Dirichlet")}
        }

    def set_initial_conditions(self, variables):
        c_inf = pybamm.Parameter("Bulk electrolyte solvent concentration [mol.m-3]")
        L_0 = pybamm.Parameter("Initial thickness [m]")
        L = variables["SEI thickness [m]"]
        c_init = c_inf
        L_init = L_0
        c = variables["Solvent concentration [mol.m-3]"]
        self.initial_conditions = {c: c_init, L: L_init}