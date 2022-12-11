import numpy as np
from casatools import table
from dataclasses import dataclass, field
from pathlib import Path
import astropy.units as un
from typing import Union
from astropy.units import Quantity
from abc import ABCMeta
from casatools import ctsys
from ..functions import FluxFunction, PolFunction

tb = table()


@dataclass(init=True, repr=True)
class PolarizedSource(metaclass=ABCMeta):
    # Object that takes information of different known
    # polarized sources from https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol
    nu: Quantity = None
    flux_density: Quantity = None
    pol_angle: Quantity = None
    pol_fraction: np.ndarray = None
    spectral_idx_coefficients: np.ndarray = None
    spectral_idx_coefficients_errors: np.ndarray = None
    source: str = None
    spix_dict: dict = field(default_factory=dict)

    def __post_init__(self):
        self.spix_dict = {
            'Perley-Butler 2013': 'PerleyButler2013Coeffs',
            'Perley-Butler 2017': 'PerleyButler2017Coeffs',
            'Scaife-Heald 2012': 'ScaifeHeald2012Coeffs'
        }
        switcher = {
            "3C48": "p3c48",
            "3C48_2019": "p3c48_2019",
            "3C138": "p3c138",
            "3C138_2019": "p3c138_2019",
            "3C286": "p3c286",
            "3C286_2019": "p3c286_2019",
            "3C147": "p3c147",
            "3C147_2019": "p3c147_2019",
            "3C295_2019": "p3c295_2019",
            "3C196_2019": "p3c196_2019"
        }

        if self.source is not None:
            self.source = self.source.upper()

        if self.source != "":
            method = getattr(self, switcher.get(self.source, "__init_empty"))
            method()

        # convert frequencies to Hz
        self.nu = self.nu.to(un.Hz)

        # Pol angle from degrees to radians
        self.pol_angle = self.pol_angle.to(un.rad)

        # Pol fraction in percentage to fraction
        self.pol_fraction /= 100.0

    def p3c48(self):
        self.nu = np.array(
            [
                1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80,
                12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5
            ]
        ) * un.GHz
        self.pol_angle = np.array(
            [
                25.0, 140.0, -5.0, -150.0, -120.0, -100.0, -92.0, -84.0, -75.0, -72.0, -68.0, -67.0,
                -64.0, -62.0, -62.0, -62.0, -63.0, -64.0, -66.0, -67.0, -70.0, -70.0, -77.0, -85.0
            ]
        ) * un.deg
        self.pol_fraction = np.array(
            [
                0.3, 0.5, 0.7, 0.9, 1.4, 2.0, 2.5, 3.2, 3.8, 4.2, 5.2, 5.2, 5.3, 5.4, 6.0, 6.1, 6.4,
                6.4, 6.9, 7.1, 7.7, 7.8, 7.4, 7.5
            ]
        )

        self.flux_density = np.zeros_like(self.pol_fraction) * un.Jy

    def p3c48_2019(self):
        self.__read_source_from_txt()
        self.source = "3C48"

    def p3c138(self):
        self.nu = np.array(
            [
                1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80,
                12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5
            ]
        ) * un.GHz
        self.pol_angle = np.array(
            [
                -14.0, -11.0, -10.0, -10.0, -9.0, -10.0, -10.0, 0.0, -11.0, -11.0, -12.0, -12.0,
                -10.0, -8.0, -7.0, -7.0, -8.0, -9.0, -12.0, -13.0, -16.0, -17.0, -24.0, -27.0
            ]
        ) * un.deg
        self.pol_fraction = np.array(
            [
                5.6, 7.5, 8.4, 9.0, 10.4, 10.7, 10.0, 0.0, 10.0, 10.4, 9.8, 10.0, 10.4, 10.1, 8.4,
                7.9, 7.7, 7.4, 6.7, 6.5, 6.7, 6.6, 6.6, 6.5
            ]
        )
        self.flux_density = np.zeros_like(self.pol_fraction) * un.Jy

    def p3c138_2019(self):
        self.__read_source_from_txt()
        self.source = "3C138"

    def p3c286(self):
        self.nu = np.array(
            [
                1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80,
                12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5
            ]
        ) * un.GHz
        self.pol_angle = np.array(
            [
                33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 33.0, 34.0, 34.0,
                34.0, 34.0, 34.0, 34.0, 34.0, 35.0, 35.0, 35.0, 36.0, 36.0
            ]
        ) * un.deg
        self.pol_fraction = np.array(
            [
                8.6, 9.5, 9.9, 10.1, 10.5, 10.8, 10.9, 11.1, 11.3, 11.4, 11.6, 11.7, 11.9, 11.9,
                11.9, 11.9, 12.1, 12.2, 12.5, 12.5, 12.6, 12.6, 13.1, 13.2
            ]
        )
        self.flux_density = np.zeros_like(self.pol_fraction) * un.Jy

    def p3c286_2019(self):
        self.__read_source_from_txt()
        self.source = "3C286"

    def p3c147(self):
        self.nu = np.array(
            [
                1.05, 1.45, 1.64, 1.95, 2.45, 2.95, 3.25, 3.75, 4.50, 5.00, 6.50, 7.25, 8.10, 8.80,
                12.8, 13.7, 14.6, 15.5, 18.1, 19.0, 22.4, 23.3, 36.5, 43.5
            ]
        ) * un.GHz
        self.pol_angle = np.array(
            [
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0, 0.0, -65.0, -39.0, -24.0, -11.0,
                43.0, 48.0, 53.0, 59.0, 67.0, 68.0, 75.0, 76.0, 85.0, 86.0
            ]
        ) * un.deg
        self.pol_fraction = np.array(
            [
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.3, 0.6, 0.7, 0.8, 2.2, 2.4, 2.7,
                2.9, 3.4, 3.5, 3.8, 3.8, 4.4, 5.2
            ]
        )
        self.flux_density = np.zeros_like(self.pol_fraction) * un.Jy

    def p3c147_2019(self):
        self.__read_source_from_txt()
        self.source = "3C147"

    def p3c196_2019(self):
        self.__read_source_from_txt()
        self.source = "3C196"

    def p3c295_2019(self):
        self.__read_source_from_txt()
        self.source = "3C295"

    def __read_source_from_txt(self):
        path_txt_files = Path(__file__).parent.parent / "sources_2019"
        file_name = path_txt_files / (self.source.lower() + ".txt")
        data = np.loadtxt(str(file_name))
        self.nu = data[:, 0] * un.GHz
        self.flux_density = data[:, 1] * un.Jy
        self.pol_fraction = data[:, 2]
        self.pol_angle = data[:, 3] * un.rad

    def __init_empty(self):
        self.nu = np.array([]) * un.GHz
        self.pol_angle = np.array([]) * un.deg
        self.pol_fraction = np.array([])

    def is_correct(self):
        if len(self.nu) == len(self.pol_angle) and len(self.nu) == len(self.pol_fraction):
            return True
        else:
            return False

    # Returns values in an array lower than a certain frequency
    @staticmethod
    def filter(nu, data, nu_min: [float, Quantity] = 0.0, nu_max: [float, Quantity] = np.inf):
        valid_indexes = np.where((nu >= nu_min) & (nu <= nu_max))
        nu = nu[valid_indexes]
        data = data[valid_indexes]
        return nu, data

    @staticmethod
    def flux_scalar_giving_coefficients(
        nu: Union[float, Quantity], coefficients: np.ndarray
    ) -> float:
        # To use this formula frequency nu must be in GHz
        if isinstance(nu, Quantity):
            nu = nu.to(un.GHz).value

        index = np.arange(0, len(coefficients))
        flux_at_nu = np.sum(coefficients * np.log10(nu)**index)

        return 10.0**flux_at_nu

    @staticmethod
    def flux_giving_coefficients(
        nu: Union[np.ndarray, Quantity], coefficients: np.ndarray
    ) -> np.ndarray:
        # To use this formula frequency nu must be in GHz
        if isinstance(nu, Quantity):
            nu = nu.to(un.GHz).value
        flux_at_nu = np.zeros_like(nu, dtype=np.float32)

        for i in range(len(coefficients)):
            flux_at_nu += coefficients[i] * np.log10(nu)**i

        return 10.0**flux_at_nu

    def flux_scalar(self, nu: Union[float, Quantity]) -> float:
        return self.flux_scalar_giving_coefficients(nu, self.spectral_idx_coefficients)

    def flux(self, nu: Union[np.ndarray, Quantity]) -> np.ndarray:
        return self.flux_giving_coefficients(nu, self.spectral_idx_coefficients)

    def get_coefficients_from_table(self, standard="Perley-Butler 2017", epoch="2017"):
        coefficients_table = ctsys.resolve("nrao/VLA/standards/") + self.spix_dict[standard]
        tb.open(coefficients_table)
        _query_table = tb.taql("select * from " + coefficients_table + " where Epoch=" + epoch)
        coefficients = _query_table.getcol(self.source + "_coeffs").flatten()
        coefficients_errs = _query_table.getcol(self.source + "_coefferrs").flatten()
        if coefficients.size == 0 or coefficients_errs.size == 0:
            raise ValueError("The selected epoch does not have any data")

        self.spectral_idx_coefficients = coefficients
        self.spectral_idx_coefficients_errors = coefficients_errs
        tb.close()
        return coefficients

    def get_coefficients_from_flux(self, nu: Quantity, nu_0: Quantity = None):
        if nu_0 is None:
            nu_0 = (np.max(nu) + np.min(nu)) / 2.

        nearest_nu_0_index = np.argmin(np.abs(nu - nu_0))
        flux_0 = self.flux_density[nearest_nu_0_index]
        initial_coefficients = np.random.rand(2)
        source_func_frac = FluxFunction(flux_0=flux_0, xdata=nu, x_0=nu_0)
        source_func_frac.fit(nu, self.flux_density.value, initial_coefficients)
        return source_func_frac.coefficients, source_func_frac.coefficients_errors

    def fit_alpha_and_beta(self, nu: Quantity, nu_0: Quantity = None):
        if nu_0 is None:
            nu_0 = (np.max(nu) + np.min(nu)) / 2.
        flux_0 = self.flux_scalar(nu_0)
        fluxes = self.flux(nu)
        upper_bound = self.spectral_idx_coefficients + self.spectral_idx_coefficients_errors
        lower_bound = self.spectral_idx_coefficients - self.spectral_idx_coefficients_errors
        fluxes_upper_bound = self.flux_giving_coefficients(nu, upper_bound)
        fluxes_lower_bound = self.flux_giving_coefficients(nu, lower_bound)
        error_sigma = 0.5 * (fluxes_upper_bound - fluxes_lower_bound)
        if np.sum(error_sigma) == 0.0:
            error_sigma = None

        initial_coefficients = np.random.rand(2)
        source_func_frac = FluxFunction(flux_0=flux_0, xdata=nu, x_0=nu_0)
        source_func_frac.fit(nu, fluxes, initial_coefficients, sigma=error_sigma)
        return source_func_frac.coefficients, source_func_frac.coefficients_errors

    # Returns pol fraction coeffs
    def get_pol_fraction_coefficients(
        self,
        n_terms: int = 3,
        nu_0: [float, Quantity] = None,
        nu_min: [float, Quantity] = 0.0,
        nu_max: [float, Quantity] = np.inf
    ):
        nu, pol_frac = self.filter(self.nu, self.pol_fraction, nu_min, nu_max)
        if nu_0 is None:
            nu_0 = (np.max(nu) + np.min(nu)) / 2.
        initial_pol_fraction_coefficients = np.random.uniform(0.0, 1.0, n_terms)
        source_func_frac = PolFunction(x_0=nu_0, n_terms=n_terms)
        source_func_frac.fit(nu, pol_frac, initial_pol_fraction_coefficients)
        return source_func_frac.coefficients, source_func_frac.coefficients_errors

    # Returns pol angle coefficients in radians
    def get_pol_angle_coefficients(
        self,
        n_terms: int = 3,
        nu_0: [float, Quantity] = None,
        nu_min: [float, Quantity] = 0.0,
        nu_max: [float, Quantity] = np.inf
    ):
        nu, pol_angle = self.filter(self.nu, self.pol_angle, nu_min, nu_max)
        if nu_0 is None:
            nu_0 = (np.max(nu) + np.min(nu)) / 2.
        initial_pol_angle_coefficients = np.random.uniform(-np.pi, np.pi, n_terms)
        source_func_angle = PolFunction(x_0=nu_0, n_terms=n_terms)
        source_func_angle.fit(nu, pol_angle, initial_pol_angle_coefficients)
        return source_func_angle.coefficients, source_func_angle.coefficients_errors

    def get_known_source_information(
        self, nu_0: Quantity = 0.0, standard: str = "Perley-Butler 2017", epoch: str = "2017"
    ):
        self.get_coefficients_from_table(standard=standard, epoch=epoch)
        nu_fit = np.linspace(0.3275, 50.0, 40) * un.GHz
        spec_idx, spec_idx_err = self.fit_alpha_and_beta(nu_fit, nu_0=nu_0)
        intensity = self.flux_scalar(nu_0)
        return intensity, spec_idx, spec_idx_err

    def get_source_polarization_information(
        self, n_terms_angle=3, n_terms_frac=3, nu_min=0.0, nu_max=np.inf
    ):
        pol_fraction_coefficients, pol_fraction_coefficients_errors = self.get_pol_fraction_coefficients(
            n_terms=n_terms_frac, nu_min=nu_min, nu_max=nu_max
        )
        pol_angle_coefficients, pol_angle_coefficients_errors = self.get_pol_angle_coefficients(
            n_terms=n_terms_angle, nu_min=nu_min, nu_max=nu_max
        )
        return pol_angle_coefficients, pol_angle_coefficients_errors, pol_fraction_coefficients, pol_fraction_coefficients_errors
