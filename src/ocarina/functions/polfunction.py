from .function import Function
from dataclasses import dataclass, field
import numpy as np


@dataclass(init=True, repr=True)
class PolFunction(Function):
    n_terms: int = 0
    coefficients: np.ndarray = None

    def __post_init__(self):
        if self.coefficients is None:
            self.coefficients = np.array([])

    def f(self, xdata, *args):
        y = np.zeros_like(xdata)
        for i in range(0, self.n_terms):
            y += args[i] * np.power((xdata - self.x_0) / self.x_0, i)
        return y

    def f_eval(self, xdata, coefficients):
        y = np.zeros(len(xdata))
        for i in range(0, len(coefficients)):
            y += coefficients[i] * np.power((xdata - self.x_0) / self.x_0, i)
        return y
