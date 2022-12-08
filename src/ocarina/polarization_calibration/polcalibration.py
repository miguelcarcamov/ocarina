import os
import sys
import numpy as np
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from abc import ABCMeta
from casatasks import polcal, applycal, gaincal, setjy, fluxscale, rmtables, flagdata, flagmanager
from casaplotms import plotms
import astropy.units as un
from astropy.units import Quantity
from ..utils import query_table
from ..polarized_sources import PolarizedSource


@dataclass(init=True, repr=True)
class PolarizationCalibrator(metaclass=ABCMeta):
    vis_name: str = ""
    spw_ids: np.ndarray = None
    antennas: str = ""
    pol_angle_field: str = ""
    leakage_field: str = ""
    target: str = ""
    ref_ant: str = ""
    k_cross_ref_ant: str = ""
    mapped_spw: int = 0
    nu_0: Quantity = None
    nu_min: Quantity = None
    nu_max: Quantity = None
    old_vla: bool = False
    number_spectral_windows: int = dataclass_field(init=False, repr=True, default=0)
    k_cross_table: str = ""
    leakage_table: str = ""
    pol_angle_table: str = ""

    def __post_init__(self):

        if self.nu_0 is None:
            spw_table = query_table(
                table_name=self.vis_name,
                query="SELECT REF_FREQUENCY FROM " + self.vis_name + "/SPECTRAL_WINDOW" +
                " WHERE !FLAG_ROW"
            )
            spw_ref_frequencies = spw_table.getcol("REF_FREQUENCY") * un.Hz
            self.nu_0 = np.median(spw_ref_frequencies)

        if self.nu_min is None:
            min_freq = query_table(
                table_name=self.vis_name,
                query="SELECT GMIN(CHAN_FREQ) AS FREQ_MIN FROM " + self.vis_name +
                "/SPECTRAL_WINDOW" + " WHERE !FLAG_ROW"
            )
            self.nu_min = min_freq.getcol("FREQ_MIN")[0] * un.Hz

        if self.nu_max is None:
            max_freq = query_table(
                table_name=self.vis_name,
                query="SELECT GMAX(CHAN_FREQ) AS FREQ_MAX FROM " + self.vis_name +
                "/SPECTRAL_WINDOW" + " WHERE !FLAG_ROW"
            )
            self.nu_max = max_freq.getcol("FREQ_MAX")[0] * un.Hz

        if self.spw_ids is None:
            spw_table = query_table(
                table_name=self.vis_name,
                query="SELECT ROWID() AS ROW_ID FROM " + self.vis_name + "/SPECTRAL_WINDOW" +
                " WHERE !FLAG_ROW"
            )
            spw_dict = spw_table.getvarcol("ROW_ID")
            self.spw_ids = np.array([item for item in spw_dict.values()]).flatten()

        self.number_spectral_windows = len(self.spw_ids)
        print("Number of spectral windows: " + str(self.number_spectral_windows))
        print("Reference freq: {0}".format(self.nu_0.to(un.GHz)))
        print("Minimum freq: {0}".format(self.nu_min.to(un.GHz)))
        print("Maximum freq: {0}".format(self.nu_max.to(un.GHz)))

    def plot_models(self, field: str = ""):
        plotms(
            vis=self.vis_name,
            field=field,
            correlation='RR',
            timerange='',
            antenna=self.ref_ant,
            xaxis='frequency',
            yaxis='amp',
            ydatacolumn='model',
            showgui=False,
            plotfile=field + '_RRamp_model.png',
            overwrite=True
        )
        plotms(
            vis=self.vis_name,
            field=field,
            correlation='RL',
            timerange='',
            antenna=self.ref_ant,
            xaxis='frequency',
            yaxis='amp',
            ydatacolumn='model',
            showgui=False,
            plotfile=field + '_RLamp_model.png',
            overwrite=True
        )
        plotms(
            vis=self.vis_name,
            field=field,
            correlation='RR',
            timerange='',
            antenna=self.ref_ant,
            xaxis='frequency',
            yaxis='phase',
            ydatacolumn='model',
            showgui=False,
            plotfile=field + '_RRphase_model.png',
            overwrite=True
        )
        plotms(
            vis=self.vis_name,
            field=field,
            correlation='RL',
            timerange='',
            antenna=self.ref_ant,
            xaxis='frequency',
            yaxis='phase',
            ydatacolumn='model',
            showgui=False,
            plotfile=field + '_RLphase_model.png',
            overwrite=True
        )

    def set_unknown_model(
        self,
        pol_source_object: PolarizedSource = None,
        field: str = "",
        gain_table: str = "",
        reference_field: str = "",
        transfer_field: str = "",
        fit_order: int = 1,
        use_scratch: bool = False
    ):

        field_table = query_table(
            table_name=self.vis_name, query="SELECT NAME FROM " + self.vis_name + "/FIELD"
        )
        field_ids = field_table.rownumbers()
        fields = field_table.getcol("NAME")
        field_id_query = np.where(fields == field)[0][0]
        field_id = field_ids[field_id_query]
        print("Field " + field + " - ID: " + str(field_id))
        field_table.close()
        # TODO: get MS string without using slicing
        flux_table = self.vis_name[:-3] + ".F." + field
        if os.path.exists(flux_table):
            rmtables(flux_table)

        # From fluxscale documentation we know that coefficients are
        # returned from the natural log nu/nu_0 Taylor expansion
        flux_dict = fluxscale(
            vis=self.vis_name,
            fluxtable=flux_table,
            caltable=gain_table,
            reference=reference_field,
            transfer=transfer_field,
            fitorder=fit_order
        )
        print(flux_dict)
        coefficients = np.array(flux_dict[str(field_id)]['spidx'].tolist())

        print(
            "Coefficients: ", coefficients
        )  # a0 log10(S at nu_0), a1 spectral idx, a2 spectral curvature
        pol_source_object.coefficients = coefficients
        # Extract a0 and make coefficients to have only spectral index and spectral curvature coefficients
        intensity = 10.0**coefficients[0]
        spectral_index = coefficients[1:-1]

        print("Setting model of: " + pol_source_object.source)
        print("Field: " + field)
        print("Reference freq: {0}".format(self.nu_0.to(un.GHz)))
        print("I(nu_0) = " + str(intensity))

        print("Alpha & Beta: ", spectral_index)

        source_dict = setjy(
            vis=self.vis_name,
            field=field,
            standard='manual',
            spw='',
            fluxdensity=[intensity, 0, 0, 0],
            spix=spectral_index.tolist(),
            reffreq=str(self.nu_0),
            interpolation="nearest",
            scalebychan=True,
            usescratch=use_scratch
        )
        print(source_dict)
        self.plot_models(field)

        return flux_table

    def set_known_model(
        self,
        pol_source_object: PolarizedSource = None,
        standard: str = "Perley-Butler 2017",
        field: str = "",
        epoch: str = "2017",
        n_terms_angle: int = 3,
        n_terms_frac: int = 3,
        use_scratch: bool = False
    ):

        # get spectral idx coeffs from NRAO VLA tables
        intensity, spec_idx, spec_idx_err = pol_source_object.get_known_source_information(
            nu_0=self.nu_0, standard=standard, epoch=epoch
        )

        pol_angle_coefficients, pol_angle_coefficients_errors, pol_fraction_coefficients, pol_fraction_coefficients_errors = pol_source_object.get_source_polarization_information(
            n_terms_angle=n_terms_angle,
            n_terms_frac=n_terms_frac,
            nu_min=self.nu_min,
            nu_max=self.nu_max
        )
        # get intensity in reference frequency
        print("Setting model of: " + pol_source_object.source)
        print("Field: " + field)
        print("Reference freq: {0}".format(self.nu_0.to(un.GHz)))
        print("I = " + str(intensity))

        print("Alpha & Beta: ", spec_idx)
        print("Error: ", spec_idx_err)
        print("Pol fraction coeffs: ", pol_fraction_coefficients)
        print("Error: ", pol_fraction_coefficients_errors)
        print("Pol angle coeffs: ", pol_angle_coefficients)
        print("Error: ", pol_angle_coefficients_errors)
        source_dict = setjy(
            vis=self.vis_name,
            field=field,
            standard='manual',
            spw='',
            fluxdensity=[intensity, 0, 0, 0],
            spix=spec_idx.tolist(),
            reffreq=str(self.nu_0),
            polindex=pol_fraction_coefficients,
            polangle=pol_angle_coefficients,
            interpolation="nearest",
            scalebychan=True,
            usescratch=use_scratch
        )
        print(source_dict)
        self.plot_models(field)

    def solve_cross_hands_delay(
        self,
        min_snr: float = 3.0,
        sol_int: str = 'inf',
        combine: str = 'scan,spw',
        spw_interval: str = "",
        channels: str = "",
        ref_ant_mode: str = "strict"
    ):
        print("Solving Cross-hand Delays")
        print("Vis: " + self.vis_name)
        print("Field: " + self.pol_angle_field)
        print("Refant: " + self.ref_ant)

        # TODO: Get string before .ms without using slicing
        cal_table = self.vis[:-3] + ".Kcross"
        if os.path.exists(cal_table):
            rmtables(cal_table)
        first_spw = self.spw_ids[0]
        last_spw = self.spw_ids[-1]

        if spw_interval == "":
            if channels == "":
                spw = str(first_spw) + '~' + str(last_spw)
            else:
                spw = str(first_spw) + '~' + str(last_spw) + ':' + channels
        else:
            if channels != "":
                spw = spw_interval + ':' + channels
            else:
                spw = spw_interval

        print("Spw: " + spw)

        gaincal(
            vis=self.vis_name,
            caltable=cal_table,
            field=self.pol_angle_field,
            spw=spw,
            refant=self.k_cross_ref_ant,
            refantmode=ref_ant_mode,
            antenna=self.antennas,
            minsnr=min_snr,
            gaintype="KCROSS",
            solint=sol_int,
            combine=combine,
            calmode="ap",
            append=False,
            gaintable=[''],
            gainfield=[''],
            interp=[''],
            spwmap=[[]],
            parang=True
        )

        if not os.path.exists(cal_table):
            raise FileNotFoundError("Caltable was not created and cannot continue. Exiting...")

        plotms(
            vis=cal_table,
            xaxis='frequency',
            yaxis='delay',
            antenna=self.k_cross_ref_ant,
            coloraxis='corr',
            showgui=False,
            plotfile=self.vis_name[:-3] + '.freqvsdelayKcross.png',
            overwrite=True
        )

        self.k_cross_table = cal_table
        return cal_table

    def calibrate_leakage(
        self,
        sol_int: str = 'inf',
        min_snr: float = 3.0,
        pol_type: str = "Df",
        spw_map: list = [],
        gain_table: list = [],
        gain_field: list = [],
        clip_min: float = 0.0,
        clip_max: float = 0.25,
        flag_clip: bool = True,
        interp_mode: str = 'linear',
        spw: str = "",
        field: str = ""
    ):
        if not gain_table:
            if self.k_cross_table == "":
                gain_table = []
            else:
                gain_table = [self.k_cross_table]
        print("Leakage calibration")
        print("Vis: " + self.vis_name)

        print("Gain tables: ", gain_table)
        print("Refant: " + self.ref_ant)

        if field == "":
            cal_table = self.vis_name[:-3] + ".D0"
            print("Field: " + self.leakage_field)
        else:
            cal_table = self.vis[:-3] + ".D." + field
            print("Field " + field)

        if not gain_field:
            gain_field = [''] * len(gain_table)
        if os.path.exists(cal_table):
            rmtables(cal_table)

        first_spw = self.spw_ids[0]
        last_spw = self.spw_ids[-1]

        if spw == "":
            spw = str(first_spw) + '~' + str(last_spw)

        spw_map0 = [self.mapped_spw] * self.number_spectral_windows

        if not spw_map:
            if len(gain_table) - 1 > 0:
                spw_map_empty = [[]] * (len(gain_table) - 1)  # subtract kcrosstable
                spw_map_empty.insert(0, spw_map0)
                spw_map = spw_map_empty
            else:
                spw_map = [spw_map0]

        interp = [interp_mode] * len(gain_table)

        if self.old_vla:
            spw = ''
            spw_map = []
            interp = 'nearest'

        print("Spw: " + spw)
        print("Spwmap: ", spw_map)

        if field == "":
            field = self.leakage_field

        polcal(
            vis=self.vis_name,
            caltable=cal_table,
            field=field,
            spw=spw,
            refant=self.ref_ant,
            antenna=self.antennas,
            poltype=pol_type,
            solint=sol_int,
            spwmap=spw_map,
            combine='scan',
            interp=interp,
            minsnr=min_snr,
            gaintable=gain_table,
            gainfield=gain_field
        )

        if not os.path.exists(cal_table):
            raise FileNotFoundError("Caltable was not created and cannot continue. Exiting...")

        if flag_clip:
            flagdata(
                vis=cal_table,
                mode='clip',
                correlation='ABS_ALL',
                clipminmax=[clip_min, clip_max],
                datacolumn='CPARAM',
                clipoutside=True,
                action='apply',
                flagbackup=False,
                savepars=False
            )
            flagmanager(
                vis=cal_table,
                mode="save",
                versionname="clip_flagging",
                comment="Clip flagging outside [" + str(clip_min) + "," + str(clip_max) + "]"
            )

        plotms(
            vis=cal_table,
            xaxis='freq',
            yaxis='amp',
            iteraxis='antenna',
            coloraxis='corr',
            showgui=False,
            plotfile=self.vis_name[:-3] + '.D0.ampvsfreq.png',
            overwrite=True
        )

        plotms(
            vis=cal_table,
            xaxis='chan',
            yaxis='phase',
            iteraxis='antenna',
            coloraxis='corr',
            plotrange=[-1, -1, -180, 180],
            showgui=False,
            plotfile=self.vis_name[:-3] + '.D0.phasevschan.png',
            overwrite=True
        )

        plotms(
            vis=cal_table,
            xaxis='chan',
            yaxis='phase',
            iteraxis='antenna',
            coloraxis='corr',
            plotrange=[-1, -1, -180, 180],
            showgui=False,
            plotfile=self.vis_name[:-3] + '.D0.ampvsantenna.png',
            overwrite=True
        )

        plotms(
            vis=cal_table,
            xaxis='antenna1',
            yaxis='amp',
            coloraxis='corr',
            showgui=False,
            plotfile=self.vis_name[:-3] + '.D0.ampvsantenna1.png',
            overwrite=True
        )

        self.leakage_table = cal_table
        return cal_table

    def calibrate_pol_angle(
        self,
        sol_int: str = 'inf',
        min_snr: float = 3.0,
        pol_type: str = "Xf",
        spw_map: list = [],
        gain_table: list = [],
        gain_field: list = [],
        interp_mode: str = 'linear',
        spw: str = "",
        field: str = ""
    ):
        if not gain_table:
            if self.k_cross_table == "":
                gain_table = [self.leakage_table]
            else:
                gain_table = [self.k_cross_table, self.leakage_table]
        print("Polarization angle calibration")
        print("Vis: " + self.vis_name)
        print("Field: " + self.pol_angle_field)

        print("Gain tables: ", gain_table)

        print("Refant: " + self.ref_ant)

        # TODO: fix string
        cal_table = self.vis_name[:-3] + ".X0"
        if not gain_field:
            gain_field = [''] * len(gain_table)

        if os.path.exists(cal_table):
            rmtables(cal_table)
        first_spw = self.spw_ids[0]
        last_spw = self.spw_ids[-1]

        if spw == "":
            spw = str(first_spw) + '~' + str(last_spw)

        spw_map0 = [self.mapped_spw] * self.number_spectral_windows

        if not spw_map:
            if len(gain_table) - 1 > 0:
                spw_map_empty = [[]] * (len(gain_table) - 1)  # subtract kcrosstable
                spw_map_empty.insert(0, spw_map0)
                spw_map = spw_map_empty
            else:
                spw_map = [spw_map0]

        interp = [interp_mode] * len(gain_table)
        if self.old_vla:
            spw = ''
            spw_map = []
            interp = 'nearest'

        print("Spw: " + spw)
        print("Spwmap: ", spw_map)

        if field == "":
            field = self.pol_angle_field

        polcal(
            vis=self.vis_name,
            caltable=cal_table,
            field=field,
            spw=spw,
            refant=self.ref_ant,
            antenna=self.antennas,
            poltype=pol_type,
            solint=sol_int,
            combine='scan',
            spwmap=spw_map,
            interp=interp,
            minsnr=min_snr,
            gaintable=gain_table,
            gainfield=gain_field
        )

        if not os.path.exists(cal_table):
            sys.exit("Caltable was not created and cannot continue. Exiting...")
        plotms(
            vis=cal_table,
            xaxis='frequency',
            yaxis='phase',
            coloraxis='spw',
            showgui=False,
            plotfile=self.vis_name[:-3] + '.X0.phasevsfreq.png',
            overwrite=True
        )

        self.pol_angle_table = cal_table
        return cal_table

    def plot_leakage(self, plot_dir="", field="", cal_table=""):
        if field == "" and cal_table == "":
            plotms(
                vis=self.leakage_table,
                xaxis='antenna',
                yaxis='amp',
                plotfile=plot_dir + self.vis_name[:-3] + '.D0.amp.png',
                showgui=False,
                overwrite=True
            )
            plotms(
                vis=self.leakage_table,
                xaxis='antenna',
                yaxis='phase',
                iteraxis='antenna',
                plotfile=plot_dir + self.vis_name[:-3] + '.D0.phs.png',
                showgui=False,
                overwrite=True
            )
            plotms(
                vis=self.leakage_table,
                xaxis='antenna',
                yaxis='snr',
                showgui=False,
                plotfile=plot_dir + self.vis_name[:-3] + '.D0.snr.png',
                overwrite=True
            )
            plotms(
                vis=self.leakage_table,
                xaxis='real',
                yaxis='imag',
                showgui=False,
                plotfile=plot_dir + self.vis_name[:-3] + '.D0.cmplx.png',
                overwrite=True
            )
        else:
            plotms(
                vis=cal_table,
                xaxis='antenna',
                yaxis='amp',
                plotfile=plot_dir + self.vis_name[:-3] + '.D.' + field + 'amp.png',
                showgui=False,
                overwrite=True
            )
            plotms(
                vis=cal_table,
                xaxis='antenna',
                yaxis='phase',
                iteraxis='antenna',
                plotfile=plot_dir + self.vis_name[:-3] + '.D.' + field + 'phs.png',
                showgui=False,
                overwrite=True
            )
            plotms(
                vis=cal_table,
                xaxis='antenna',
                yaxis='snr',
                showgui=False,
                plotfile=plot_dir + self.vis_name[:-3] + '.D.' + field + 'snr.png',
                overwrite=True
            )
            plotms(
                vis=cal_table,
                xaxis='real',
                yaxis='imag',
                showgui=False,
                plotfile=plot_dir + self.vis_name[:-3] + '.D.' + field + 'cmplx.png',
                overwrite=True
            )

    def apply_single_solution(
        self,
        field: str = '',
        spw: str = '',
        gain_table: list = [],
        gain_field: list = [],
        select_data: bool = True,
        spw_map: list = [],
        cal_wt: list = [False],
        apply_mode: str = "calflagstrict",
        interp: str = 'linear',
        antenna: str = '',
        flag_backup: bool = True
    ):
        applycal(
            vis=self.vis_name,
            field=field,
            spw=spw,
            gaintable=gain_table,
            gainfield=gain_field,
            selectdata=select_data,
            spwmap=spw_map,
            calwt=cal_wt,
            applymode=apply_mode,
            interp=interp,
            antenna=antenna,
            parang=True,
            flagbackup=flag_backup
        )

    def apply_solutions(
        self,
        spw_map: list = [],
        gain_table: list = [],
        gain_field: list = [],
        apply_mode: str = "calflagstrict",
        antenna: str = '*&*',
        flag_backup: bool = True
    ):
        if not gain_table:
            if self.k_cross_table == "":
                gain_table = [self.leakage_table, self.pol_angle_table]
            else:
                gain_table = [self.k_cross_table, self.leakage_table, self.pol_angle_table]
        print("Applying solutions")
        print("Gain tables: ", gain_table)
        first_spw = self.spw_ids[0]
        last_spw = self.spw_ids[-1]

        interp = [''] * len(gain_table)
        spw = str(first_spw) + '~' + str(last_spw)

        spw_map0 = [self.mapped_spw] * self.number_spectral_windows

        if not spw_map:
            if len(gain_table) - 1 > 0:
                spw_map_empty = [[]] * (len(gain_table) - 1)  # subtract kcrosstable
                spw_map_empty.insert(0, spw_map0)
                spw_map = spw_map_empty
            else:
                spw_map = [spw_map0]

        cal_wt = [False] * len(gain_table)
        select_data = True
        if self.old_vla:
            interp = 'nearest'
            spw_map = []
            spw = ''
            cal_wt = [False]
            select_data = False
            antenna = ''

        print("Spw: " + spw)
        print("Spwmap: ", spw_map)
        if not gain_field:
            gain_field = [''] * len(gain_table)
        applycal(
            vis=self.vis_name,
            field='',
            spw=spw,
            gaintable=gain_table,
            selectdata=select_data,
            spwmap=spw_map,
            calwt=cal_wt,
            applymode=apply_mode,
            interp=interp,
            gainfield=gain_field,
            antenna=antenna,
            parang=True,
            flagbackup=flag_backup
        )

    def final_plots(self):
        plotms(
            vis=self.vis_name,
            field=self.pol_angle_field,
            correlation='',
            timerange='',
            antenna='',
            avgtime='60',
            xaxis='frequency',
            yaxis='amp',
            ydatacolumn='corrected',
            coloraxis='corr',
            plotfile=self.pol_angle_field + '.corrected-amp.png',
            showgui=False,
            overwrite=True
        )

        plotms(
            vis=self.vis_name,
            field=self.pol_angle_field,
            correlation='',
            timerange='',
            antenna='',
            avgtime='60',
            xaxis='frequency',
            yaxis='phase',
            ydatacolumn='corrected',
            plotrange=[-1, -1, -180, 180],
            coloraxis='corr',
            plotfile=self.pol_angle_field + '.corrected-phase.png',
            showgui=False,
            overwrite=True
        )

        plotms(
            vis=self.vis_name,
            field=self.leakage_field,
            correlation='',
            timerange='',
            antenna='',
            avgtime='60',
            xaxis='frequency',
            yaxis='amp',
            ydatacolumn='corrected',
            coloraxis='corr',
            plotfile=self.leakage_field + '.corrected-amp.png',
            showgui=False,
            overwrite=True
        )

        plotms(
            vis=self.vis_name,
            field=self.leakage_field,
            correlation='RR,LL',
            timerange='',
            antenna='',
            avgtime='60',
            xaxis='frequency',
            yaxis='phase',
            ydatacolumn='corrected',
            plotrange=[-1, -1, -180, 180],
            coloraxis='corr',
            plotfile=self.leakage_field + '.corrected-phase.png',
            showgui=False,
            overwrite=True
        )
