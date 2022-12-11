from .function import Function
from astropy.units import Quantity
from dataclasses import dataclass
import numpy as np


@dataclass(init=True, repr=True)
class PolFunction(Function):
    n_terms: int = 0
    coefficients: np.ndarray = None

    def __post_init__(self):
        if self.coefficients is None:
            self.coefficients = np.array([])

    def f(self, xdata, *args):
        self.check_same_units(xdata)
        y = np.zeros(xdata.shape)
        nu_div = (xdata - self.x_0) / self.x_0
        if isinstance(nu_div, Quantity):
            nu_div = nu_div.value
        for i in range(0, self.n_terms):
            y += args[i] * nu_div**i
        return y

    def f_eval(self, xdata, coefficients):
        self.check_same_units(xdata)
        y = np.zeros(xdata.shape)
        nu_div = ((xdata - self.x_0) / self.x_0)
        if isinstance(nu_div, Quantity):
            nu_div = nu_div.value
        for i in range(0, len(coefficients)):
            y += coefficients[i] * nu_div**i
        return y
