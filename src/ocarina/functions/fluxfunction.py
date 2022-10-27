import numpy as np
from dataclasses import dataclass
from .function import Function


@dataclass(init=True, repr=True)
class FluxFunction(Function):
    flux_0: float = 0.0

    def f(self, xdata, *args):
        nu_div = xdata / self.x_0
        exp = args[0] + args[1] * np.log(nu_div)
        s_flux = self.flux_0 * nu_div**exp
        return s_flux

    def f_eval(self, xdata, coefficients):
        nu_div = xdata / self.x_0
        exp = coefficients[0] + coefficients[1] * np.log(nu_div)
        s_flux = self.flux_0 * nu_div**exp
        return s_flux
