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
    def f_eval(self, xdata: Union[float, Quantity], coefficients):
        pass

    def check_same_units(self, xdata: Union[float, Quantity]):
        if isinstance(xdata, Quantity) and isinstance(self.x_0, Quantity):
            if xdata.unit != self.x_0.unit:
                print("Converting x_0 to {0}".format(xdata.unit))
                self.x_0 = self.x_0.to(xdata.unit)

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

        self.check_same_units(xdata)
        popt, pcov = curve_fit(
            self.f,
            xdata,
            data,
            p0=initial_coefficients,
            check_finite=True,
            sigma=sigma,
            bounds=(lower_bounds, upper_bounds)
        )
        opt_error = np.sqrt(np.diag(pcov))
        self.coefficients = popt
        self.coefficients_errors = opt_error
        return popt, opt_error
