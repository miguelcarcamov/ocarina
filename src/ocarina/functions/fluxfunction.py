import numpy as np
from typing import Union
from astropy.units import Quantity
from dataclasses import dataclass
from .function import Function


@dataclass(init=True, repr=True)
class FluxFunction(Function):

    def f(self, xdata: Union[np.ndarray, Quantity], *args):
        nu_div = xdata / self.x_0
        if isinstance(nu_div, Quantity):
            nu_div = nu_div.value
        exp = args[1] + args[2] * np.log10(nu_div)
        s_flux = args[0] * nu_div**exp
        return s_flux

    def f_eval(self, xdata: Union[np.ndarray, Quantity], coefficients):
        self.check_same_units(xdata)
        nu_div = xdata / self.x_0
        if isinstance(nu_div, Quantity):
            nu_div = nu_div.value
        exp = coefficients[1] + coefficients[2] * np.log10(nu_div)
        s_flux = coefficients[0] * nu_div**exp
        return s_flux
