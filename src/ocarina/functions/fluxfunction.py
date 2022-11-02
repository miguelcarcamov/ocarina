import numpy as np
from typing import Union
from astropy.units import Quantity
from dataclasses import dataclass
from .function import Function


@dataclass(init=True, repr=True)
class FluxFunction(Function):
    flux_0: float = 0.0

    def f(self, xdata: Union[np.ndarray, Quantity], *args):
        self.check_same_units(xdata)
        nu_div = xdata / self.x_0
        exp = args[0] + args[1] * np.log(nu_div)
        s_flux = self.flux_0 * nu_div**exp
        return s_flux

    def f_eval(self, xdata: Union[np.ndarray, Quantity], coefficients):
        self.check_same_units(xdata)
        nu_div = xdata / self.x_0
        exp = coefficients[0] + coefficients[1] * np.log(nu_div)
        s_flux = self.flux_0 * nu_div**exp
        return s_flux
