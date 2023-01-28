from pylab import *
import h5py
import lal
import lalsimulation as lalsim
import scipy.signal
import os
os.environ["LAL_DATA_PATH"] = os.path.join(os.environ['HOME'], "lscsoft/src/lalsuite-extra/data/lalsimulation/")
import ringdown as rd

def get_tgps_aps(tgps_geocent, ra, dec, psi=0, ifos=('H1', 'L1')):
    """ Returns times of arrival and antenna patterns at different detectors 
    for a given geocenter triggertime and source location/orientation.
    
    Arguments
    ---------
    tgps_geocent: float
        GPS time of arrival at geocenter
    ra: float
        source right ascension
    dec: float
        source declination
    psi: float
        polarization angle (def: 0)
    ifos: list of strs
        detector names (def: ['H1', 'L1'])
    
    Returns
    -------
    tgps_dict: dict
        dictionary of arrival GPS times at each detector
    ap_dict: dict
        dictionary of antenna patterns at each detector
    """
    lal_t_geo = lal.LIGOTimeGPS(tgps_geocent)
    gmst = lal.GreenwichMeanSiderealTime(lal_t_geo)
    tgps_dict = {}
    ap_dict = {}
    for ifo in ifos:
        d = lal.cached_detector_by_prefix[ifo]
        dt_ifo = lal.TimeDelayFromEarthCenter(d.location, ra, dec, lal_t_geo)
        tgps_dict[ifo] = tgps_geocent + dt_ifo
        ap_dict[ifo] = lal.ComputeDetAMResponse(d.response, ra, dec, psi, gmst)
    return tgps_dict, ap_dict

def generate_lal_hphc(approximant_key, m1_msun=None, m2_msun=None, chi1=None,
                      chi2=None, dist_mpc=None, dt=None, f_low=20, f_ref=None,
                      inclination=None, phi_ref=None, ell_max=None,
                      single_mode=None, epoch=None, mtot_msun=None, 
                      nr_path=None):

    approximant = lalsim.SimInspiralGetApproximantFromString(approximant_key)

    param_dict = lal.CreateDict()
    
    if f_ref is None:
        f_ref = f_low

    # NR handling based on https://arxiv.org/abs/1703.01076
    if approximant_key == 'NR_hdf5':
        # get masses
        mtot_msun = mtot_msun or m1_msun + m2_msun
        with h5py.File(nr_path, 'r') as f:
            m1 = f.attrs['mass1']
            m2 = f.attrs['mass2']
            m1_msun = m1 * mtot_msun/(m1 + m2)
            m2_msun = m2 * mtot_msun/(m1 + m2)
        # Compute spins in the LAL frame
        s1x, s1y, s1z, s2x, s2y, s2z = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(f_ref, mtot_msun, nr_path)
        chi1 = [s1x, s1y, s1z]
        chi2 = [s2x, s2y, s2z]
        # Create a dictionary and pass /PATH/TO/H5File
        lalsim.SimInspiralWaveformParamsInsertNumRelData(param_dict, nr_path)
        longAscNodes = np.pi / 2
    else:
        longAscNodes = 0.

    m1_kg = m1_msun*lal.MSUN_SI
    m2_kg = m2_msun*lal.MSUN_SI
    
    distance = dist_mpc*1e6*lal.PC_SI

    if single_mode is not None and ell_max is not None:
        raise Exception("Specify only one of single_mode or ell_max")

    if ell_max is not None:
        # If ell_max, load all modes with ell <= ell_max
        ma = lalsim.SimInspiralCreateModeArray()
        for ell in range(2, ell_max+1):
            lalsim.SimInspiralModeArrayActivateAllModesAtL(ma, ell)
        lalsim.SimInspiralWaveformParamsInsertModeArray(param_dict, ma)
    elif single_mode is not None:
        # If a single_mode is given, load only that mode (l,m) and (l,-m)
        param_dict = set_single_mode(param_dict, single_mode[0], single_mode[1])

    hp, hc = lalsim.SimInspiralChooseTDWaveform(m1_kg, m2_kg,
                                                chi1[0], chi1[1], chi1[2],
                                                chi2[0], chi2[1], chi2[2],
                                                distance, inclination,
                                                phi_ref, longAscNodes,
                                                0., 0., dt, f_low, f_ref,
                                                param_dict, approximant)
    return hp, hc

def generate_lal_waveform(*args, **kwargs):
    times = kwargs.pop('times')
    triggertime = kwargs.pop('triggertime')
    manual_epoch = kwargs.pop('manual_epoch', False)
    
    bufLength = len(times)
    delta_t = times[1] - times[0]
    tStart = times[0]
    tEnd = tStart + delta_t * bufLength

    kwargs['dt'] = delta_t

    hplus = kwargs.pop('hplus', None)
    hcross = kwargs.pop('hcross', None)
    if (hplus is None) or (hcross is None):
        hplus, hcross = generate_lal_hphc(*args, **kwargs)
    
    # align waveform, based on LALInferenceTemplate
    # https://git.ligo.org/lscsoft/lalsuite/blob/master/lalinference/lib/LALInferenceTemplate.c#L1124

    # /* The nearest sample in model buffer to the desired tc. */
    tcSample = round((triggertime - tStart)/delta_t)

    # /* The actual coalescence time that corresponds to the buffer
    #    sample on which the waveform's tC lands. */
    # i.e. the nearest time in the buffer
    injTc = tStart + tcSample*delta_t

    # /* The sample at which the waveform reaches tc. */
    if manual_epoch:
        # manually find peak of the waveform envelope
        habs = np.sqrt(hplus.data.data**2 + hcross.data.data**2)
        waveTcSample = np.argmax(habs)
    else:
        hplus_epoch = hplus.epoch.gpsSeconds + hplus.epoch.gpsNanoSeconds*1E-9
        waveTcSample = round(-hplus_epoch/delta_t)

    # /* 1 + (number of samples post-tc in waveform) */
    wavePostTc = hplus.data.length - waveTcSample

    # bufStartIndex = (tcSample >= waveTcSample ? tcSample - waveTcSample : 0);
    bufStartIndex = int(tcSample - waveTcSample if tcSample >= waveTcSample else 0)
    # size_t bufEndIndex = (wavePostTc + tcSample <= bufLength ? wavePostTc + tcSample : bufLength);
    bufEndIndex = int(tcSample + wavePostTc if tcSample + wavePostTc <= bufLength else bufLength)
    bufWaveLength = bufEndIndex - bufStartIndex
    waveStartIndex = int(0 if tcSample >= waveTcSample else waveTcSample - tcSample)

    if kwargs.get('window', True) and tcSample >= waveTcSample:
        # smoothly turn on waveform
        window = scipy.signal.tukey(bufWaveLength)
        window[int(0.5*bufWaveLength):] = 1.
    else:
        window = 1
    h_td = np.zeros(bufLength, dtype=complex)
    h_td[bufStartIndex:bufEndIndex] = window*hplus.data.data[waveStartIndex:waveStartIndex+bufWaveLength] -\
                                      1j*window*hcross.data.data[waveStartIndex:waveStartIndex+bufWaveLength]
    return h_td

NR_PATH = '/Users/maxisi/lscsoft/src/lvcnr-lfs/SXS/SXS_BBH_0305_Res6.h5'
def get_signal_dict(time_dict, **kwargs):
    # get antenna patterns and trigger times
    tgps_dict = kwargs.pop('tgps_dict', None)
    ap_dict = kwargs.pop('ap_dict', None)
    if not tgps_dict or not ap_dict:
        ra, dec, psi = [kwargs.pop(k) for k in ['ra', 'dec', 'psi']]
        ifos = time_dict.keys()
        tgeo = kwargs.pop('tgps_geocent')
        tgps_dict, ap_dict = get_tgps_aps(tgeo, ra, dec, psi, ifos)
        
    # get complex strain
    time = list(time_dict.values())[0]
    delta_t = time[1] - time[0]
    approx = kwargs.pop('approx')
    hp, hc = generate_lal_hphc(approx,
                               dt=delta_t, **kwargs)


    # project signal onto detector
    raw_signal_dict = {}
    for ifo, time in time_dict.items():
        h = generate_lal_waveform(hplus=hp, hcross=hc, times=time,
                                  triggertime=tgps_dict[ifo],
                                  manual_epoch=kwargs.get('manual_epoch'))
        Fp, Fc = ap_dict[ifo]
        h_ifo = Fp*h.real - Fc*h.imag

        raw_signal_dict[ifo] = h_ifo
    return raw_signal_dict, tgps_dict, ap_dict


def mcq_to_m1m2(mc,q):
    """
    Utility function for converting mchirp,q to component masses. The
    masses are defined so that m1>m2. The rvalue is a tuple (m1,m2).
    """
    factor = mc * np.power(1+q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, 2.0/5.0);
    return m1, m2

def change_spin_convention(theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2, m1, m2, f_ref, phi_orb=0.):
    iota, S1x, S1y, S1z, S2x, S2y, S2z = lalsim.SimInspiralTransformPrecessingNewInitialConditions(\
                                             theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2, \
                                             m1*lal.MSUN_SI, m2*lal.MSUN_SI, f_ref, phi_orb)

    return [S1x, S1y, S1z], [S2x, S2y, S2z], iota

def get_peak_times(*args, **kwargs):
    p = kwargs.pop('parameters')
    times = kwargs.pop('times')
    ifos = kwargs.pop('ifos', ['H1', 'L1', 'V1'])
    approx = kwargs.pop('approx', 'NRSur7dq4')

    delta_t = times[1] - times[0]
    tlen = len(times)

    fp = {k: kwargs[k] if k in kwargs else p[k] for k in ['f_ref', 'flow', 'lal_amporder']}

    chi1, chi2, iota = change_spin_convention(p['theta_jn'], p['phi_jl'], p['tilt_1'], p['tilt_2'],
                                              p['phi_12'], p['a_1'], p['a_2'], p['mass_1'], p['mass_2'],
                                              fp['f_ref'], p['phase'])

    f_start = fp['flow']*2/(fp['lal_amporder'] + 2.)
    # get strain
    h_td = generate_lal_waveform(approx, p['mass_1'], p['mass_2'], chi1, chi2, dist_mpc=p['luminosity_distance'],
                                  dt=delta_t, f_low=f_start, f_ref=fp['f_ref'], inclination=iota,
                                  phi_ref=p['phase'], ell_max=None, times=times, triggertime=p['geocent_time'])
    # FFT
    hp_td = h_td.real
    hc_td = -h_td.imag

    fft_norm = delta_t
    hp_fd = np.fft.rfft(hp_td) * fft_norm
    hc_fd = np.fft.rfft(hc_td) * fft_norm
    frequencies = np.fft.rfftfreq(tlen) / fft_norm

    # get peak time
    tp_geo_loc = np.argmax(np.abs(h_td))
    tp_geo = times[tp_geo_loc]

    geo_gps_time = lal.LIGOTimeGPS(p['geocent_time'])
    gmst = lal.GreenwichMeanSiderealTime(geo_gps_time)

    tp_dict = {'geo': tp_geo}
    for ifo in ifos:
        detector = lal.cached_detector_by_prefix[ifo]
        # get antenna patterns
        Fp, Fc = lal.ComputeDetAMResponse(detector.response, p['ra'], p['dec'], p['psi'], gmst)
        # get time delay and align waveform
        # assume reference time corresponds to envelope peak
        timedelay = lal.TimeDelayFromEarthCenter(detector.location,  p['ra'], p['dec'], geo_gps_time)

        fancy_timedelay = lal.LIGOTimeGPS(timedelay)
        timeshift = fancy_timedelay.gpsSeconds + 1e-9*fancy_timedelay.gpsNanoSeconds

        timeshift_vector = np.exp(-2.*1j*np.pi*timeshift*frequencies)

        tp_dict[ifo] = tp_geo + timedelay
    return tp_dict