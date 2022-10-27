import numpy as np
from scipy.optimize import curve_fit
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from astropy.units import Quantity
from typing import Union


@dataclass(init=True, repr=True)
class Function(metaclass=ABCMeta):

    x_0: Union[float, Quantity] = 0.0
    xdata: Union[np.ndarray, Quantity] = None
    coefficients: np.ndarray = None
    coefficients_errors: np.ndarray = None

    @abstractmethod
    def f(self, xdata, *args):
        pass

    @abstractmethod
    def f_eval(self, xdata, coefficients):
        pass

    def fit(
        self,
        xdata,
        data,
        initial_coefficients=None,
        lower_bound=-np.inf,
        upper_bound=np.inf,
        sigma=None
    ):
        lower_bounds = np.ones_like(initial_coefficients) * lower_bound
        upper_bounds = np.ones_like(initial_coefficients) * upper_bound
        popt, pcov = curve_fit(
            self.f,
            xdata,
            data,
            p0=initial_coefficients,
            check_finite=True,
            sigma=sigma,
            bounds=(lower_bounds, upper_bounds)
        )
        perr = np.sqrt(np.diag(pcov))
        self.coefficients = popt
        self.coefficients_errors = perr
        return popt, perr
