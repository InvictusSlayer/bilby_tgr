import lalsimulation
from bilby.gw.source import (_base_lal_cbc_fd_waveform,
                             _base_waveform_frequency_sequence)
from lal import CreateDict


def lal_binary_black_hole(
    frequency_array,
    mass_1,
    mass_2,
    luminosity_distance,
    a_1,
    tilt_1,
    phi_12,
    a_2,
    tilt_2,
    phi_jl,
    theta_jn,
    phase,
    dQuadMon1,
    dQuadMon2,
    **kwargs,
):
    """Introducing spin-induced quadrupole moment terms (dQuadMon1 and dQuadMon2)
        to the lal_binary_black_hole source  model using lalsimulation

    Parameters
    ----------
    frequency_array: array_like
        The frequencies at which we want to calculate the strain
    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    luminosity_distance: float
        The luminosity distance in megaparsec
    a_1: float
        Dimensionless primary spin magnitude
    tilt_1: float
        Primary tilt angle
    phi_12: float
        Azimuthal angle between the two component spins
    a_2: float
        Dimensionless secondary spin magnitude
    tilt_2: float
        Secondary tilt angle
    phi_jl: float
        Azimuthal angle between the total binary angular momentum and the
        orbital angular momentum
    theta_jn: float
        Angle between the total binary angular momentum and the line of sight
    phase: float
        The phase at coalescence
    dQuadMon1: float
        The spin-induced quadrupole moment parameter of heavier object
    dQuadMon2: float
        The spin-induced quadrupole moment parameter of lighter object
    kwargs: dict
        Optional keyword arguments
        Supported arguments:
            waveform_approximant
            reference_frequency
            minimum_frequency
            maximum_frequency
            catch_waveform_errors
            pn_spin_order
            pn_tidal_order
            pn_phase_order
            pn_amplitude_order
            mode_array:
                Activate a specific mode array and evaluate the model using those
                modes only.  e.g., waveform_arguments =
                dict(waveform_approximant='IMRPhenomHM', mode_array=[[2,2],[2,-2])
                returns the 22 and 2-2 modes only of IMRPhenomHM.  You can only
                specify modes that are included in that particular model.  e.g.,
                waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
                mode_array=[[2,2],[2,-2],[5,5],[5,-5]]) is not allowed because the
                55 modes are not included in this model.  Be aware that some models
                only take positive modes and return the positive and the negative
                mode together, while others need to call both.  e.g.,
                waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
                mode_array=[[2,2],[4,-4]]) returns the 22 and 2-2 of IMRPhenomHM.
                However, waveform_arguments =
                dict(waveform_approximant='IMRPhenomXHM', mode_array=[[2,2],[4,-4]])
                returns the 22 and 4-4 of IMRPhenomXHM.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_kwargs = dict(
        waveform_approximant="IMRPhenomPv2",
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0,
    )
    waveform_kwargs.update(kwargs)
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())
    lalsimulation.SimInspiralWaveformParamsInsertdQuadMon1(wf_dict, float(dQuadMon1))
    lalsimulation.SimInspiralWaveformParamsInsertdQuadMon2(wf_dict, float(dQuadMon2))

    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    return _base_lal_cbc_fd_waveform(
        frequency_array,
        mass_1=mass_1,
        mass_2=mass_2,
        luminosity_distance=luminosity_distance,
        a_1=a_1,
        tilt_1=tilt_1,
        phi_12=phi_12,
        a_2=a_2,
        tilt_2=tilt_2,
        phi_jl=phi_jl,
        theta_jn=theta_jn,
        phase=phase,
        dQuadMon1=dQuadMon1,
        dQuadMon2=dQuadMon2,
        **waveform_kwargs,
    )


def lal_binary_black_hole_relative_binning(
    frequency_array,
    mass_1,
    mass_2,
    luminosity_distance,
    a_1,
    tilt_1,
    phi_12,
    a_2,
    tilt_2,
    phi_jl,
    theta_jn,
    phase,
    dQuadMon1,
    dQuadMon2,
    **kwargs
):
    """ Source model to go with RelativeBinningGravitationalWaveTransient likelihood.

    All parameters are the same as in the usual source models, except `fiducial`

    fiducial: float
        If fiducial=1, waveform evaluated on the full frequency grid is returned.
        If fiducial=0, waveform evaluated at waveform_kwargs["frequency_bin_edges"]
        is returned.
    """
    fiducial = kwargs.pop("fiducial", 0)

    waveform_kwargs = dict(
        waveform_approximant='IMRPhenomPv2',
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0
    )
    waveform_kwargs.update(kwargs)
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())
    lalsimulation.SimInspiralWaveformParamsInsertdQuadMon1(wf_dict, float(dQuadMon1))
    lalsimulation.SimInspiralWaveformParamsInsertdQuadMon2(wf_dict, float(dQuadMon2))

    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    if fiducial == 1:
        _ = waveform_kwargs.pop("frequency_bin_edges", None)
        return _base_lal_cbc_fd_waveform(
            frequency_array=frequency_array,
            mass_1=mass_1,
            mass_2=mass_2,
            luminosity_distance=luminosity_distance,
            theta_jn=theta_jn,
            phase=phase,
            a_1=a_1,
            a_2=a_2,
            tilt_1=tilt_1,
            tilt_2=tilt_2,
            phi_jl=phi_jl,
            phi_12=phi_12,
            lambda_1=0.0,
            lambda_2=0.0,
            dQuadMon1=dQuadMon1,
            dQuadMon2=dQuadMon2,
            **waveform_kwargs
        )

    else:
        _ = waveform_kwargs.pop("minimum_frequency", None)
        _ = waveform_kwargs.pop("maximum_frequency", None)
        waveform_kwargs["frequencies"] = waveform_kwargs.pop("frequency_bin_edges")
        return _base_waveform_frequency_sequence(
            frequency_array=frequency_array,
            mass_1=mass_1,
            mass_2=mass_2,
            luminosity_distance=luminosity_distance,
            theta_jn=theta_jn,
            phase=phase,
            a_1=a_1,
            a_2=a_2,
            tilt_1=tilt_1,
            tilt_2=tilt_2,
            phi_jl=phi_jl,
            phi_12=phi_12,
            lambda_1=0.0,
            lambda_2=0.0,
            dQuadMon1=dQuadMon1,
            dQuadMon2=dQuadMon2,
            **waveform_kwargs
        )
