s = """# This is the input file of FORTRAN77 program "qseis06" for calculation of
# synthetic seismograms based on a layered halfspace earth model.
#
# by
# Rongjiang  Wang <wang@gfz-potsdam.de>
# GeoForschungsZentrum Potsdam
# Telegrafenberg, D-14473 Potsdam, Germany
#
# Last modified: Potsdam, Nov., 2006
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# If not specified, SI Unit System is used overall!
#
# Coordinate systems:
# cylindrical (z,r,t) with z = downward,
#                          r = from source outward,
#                          t = azmuth angle from north to east;
# cartesian (x,y,z) with   x = north,
#                          y = east,
#                          z = downward;
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#	SOURCE PARAMETERS
#	=================
# 1. source depth [km]
#------------------------------------------------------------------------------
10.00
#------------------------------------------------------------------------------
#
#	RECEIVER PARAMETERS
#	===================
# 1. receiver depth [km]
# 2. switch for distance sampling role (1/0 = equidistant/irregular); switch
#    for unit used (1/0 = km/deg)
# 3. number of distance samples
# 4. if equidistant, then start and end trace distance (> 0); else distance
#    list (please order the receiver distances from small to large)
# 5. (reduced) time begin [sec] & length of time window [sec], number of time
#    samples (<= 2*nfmax in qsglobal.h)
# 6. switch for unit of the following time reduction parameter: 1 = velocity
#    [km/sec], 0 = slowness [sec/deg]; time reduction parameter
#------------------------------------------------------------------------------
1.00
1 1
50
1.999686 99.984304
0.000000 255.500000 512
1 0.000000
#------------------------------------------------------------------------------
#
#	WAVENUMBER INTEGRATION PARAMETERS
#	=================================
# 1. select slowness integration algorithm (0 = suggested for full wave-field
#    modelling; 1 or 2 = suggested when using a slowness window with narrow
#    taper range - a technique for suppressing space-domain aliasing);
# 2. The first parameter is used to automatically estimate the wavenumber
#    truncation limit when the full-wavefield option is selected (e.g., 1e-6).
#    Smaller values yield a higher wavenumber cutoff, improving accuracy at 
#    the expense of greater computational cost. The second parameter sets the 
#    ratio of the source disk radius to the minimum epicentral distance 
#    (e.g., 0.05); larger values converge faster but deviate more from a 
#    point source.
# 3. 4 parameters for low and high slowness (Note 1) cut-offs [s/km] with
#    tapering: 0 < slw1 < slw2 defining cosine taper at the lower end, and 0 <
#    slw3 < slw4 defining the cosine taper at the higher end. default values
#    will be used in case of inconsistent input of the cut-offs (possibly with
#    much more computational effort);
# 4. parameter for sampling rate of the wavenumber integration (1 = sampled
#    with the spatial Nyquist frequency, 2 = sampled with twice higher than
#    the Nyquist, and so on: the larger this parameter, the smaller the space-
#    domain aliasing effect, but also the more computation effort); 
# 5. the factor for suppressing time domain aliasing (> 0 and <= 1) (Note 2).
#------------------------------------------------------------------------------
0
1e-6 0.05
0.000000 0.000000 0.000000 0.000000
12.000000
0.010000
#------------------------------------------------------------------------------
#
#	        OPTIONS FOR PARTIAL SOLUTIONS
#       (only applied to the source-site structure)
#	    ===========================================
#
# 1. switch for filtering free surface effects (0 = with free surface, i.e.,
#    do not select this filter; 1 = without free surface; 2 = without free
#    surface but with correction on amplitude and wave form. Note switch 2
#    can only be used for receivers at the surface)
# 2. switch for filtering waves with a shallow penetration depth (concerning
#    their whole trace from source to receiver), penetration depth limit [km]
#
#    if this option is selected, waves whose travel path never exceeds the
#    given depth limit will be filtered ("seismic nuting"). the condition for
#    selecting this filter is that the given shallow path depth limit should
#    be larger than both source and receiver depth.
#
# 3. number of depth ranges where the following selected up/down-sp2oing P or
#    SV waves should be filtered
# 4. the 1. depth range: upper and lower depth [km], switch for filtering P
#    or SV wave in this depth range:
#
#    switch no:              1      2        3       4         other
#    filtered phase:         P(up)  P(down)  SV(up)  SV(down)  Error
#
# 5. the 2. ...
#
#    The partial solution options are useful tools to increase the numerical
#    significance of desired wave phases. Especially when the desired phases
#    are smaller than the undesired phases, these options should be selected
#    and carefully combined.
#------------------------------------------------------------------------------
0
0   0.00
0
#0.0 0.0 2
#0.0 0.4 4
#------------------------------------------------------------------------------
#
#	SOURCE TIME FUNCTION (WAVELET) PARAMETERS (Note 3)
#	==================================================
# 1. wavelet duration [unit = time sample rather than sec!], that is about
#    equal to the half-amplitude cut-off period of the wavelet (> 0. if <= 0,
#    then default value = 2 time samples will be used), and switch for the
#    wavelet form (0 = user's own wavelet; 1 = default wavelet: normalized
#    square half-sinusoid for simulating a physical delta impulse; 2 = tapered
#    Heaviside wavelet, i.e. integral of wavelet 1)
# 2. IF user's own wavelet is selected, then number of the wavelet time samples
#    (<= 1024), and followed by
# 3. equidistant wavelet time samples
# 4  ...(continue) (! no comment lines allowed between the time sample list!)
#    IF default, delete line 2, 3, 4 ... or comment them out!
#------------------------------------------------------------------------------
5 2
#------------------------------------------------------------------------------
#
#	 FILTER PARAMETERS OF RECEIVERS (SEISMOMETERS OR HYDROPHONES)
#	 ============================================================
# 1. constant coefficient (normalization factor)
# 2. number of roots (<= nrootmax in qsglobal.h)
# 3. list of the root positions in the complex format (Re,Im). If no roots,
#    comment out this line
# 4. number of poles (<= npolemax in qsglobal.h)
# 5. list of the pole positions in the complex format (Re,Im). If no poles,
#    comment out this line
#------------------------------------------------------------------------------
(1.0, 0.0)
0
# (0.0, 0.0), (0.0, 0.0)
0
# (-4.35425, 4.44222), (-4.35425,-4.44222)
#------------------------------------------------------------------------------
#
#	OUTPUT FILES FOR GREEN'S FUNCTIONS (Note 4)
#	===========================================
# 1. selections of source types (yes/no = 1/0)
# 2. file names of Green's functions (please give the names without extensions,
#    which will be appended by the program automatically)
# 3. select output observables (1/0 = yes/no)
#    Note: If the wavelet form in SOURCE TIME FUNCTION (WAVELET) PARAMETERS
#    is set to a normalized square half-sinusoid, the outputs represent 
#    | velocity | volume_rate | strain_rate | stress_rate | rotation_rate | 
#    If it is set to a tapered Heaviside wavelet, the outputs represent 
#    | displacement | volume | strain | stress | rotation |
#    If a custom STF is used, the outputs are the convolution of that STF with 
#    the corresponding rate kernels.
#    File extensions of outputs:
#    *.tz, *.tr, *.tt -- velocity/disp
#    *.tv -- volume_rate/volume change (for hydrophones)
#    *.ezz, *.ezr, *.ezt, *.ett, *.err, *.ert -- stress_rate/stress
#    *.szz, *.szr, *.szt, *.stt, *.srr, *.ert -- strain_rate/strain
#    *.oz, *.or, *.ot -- rotation_rate/rotation
#------------------------------------------------------------------------------
# | explosion | strike-slip | dip-slip | clvd | single_f_v | single_f_h |
#------------------------------------------------------------------------------
1    1    1    1    0    0
'ex' 'ss' 'ds' 'cl' 'fz' 'fh'
#------------------------------------------------------------------------------
# | velocity | volume_rate | strain_rate | stress_rate | rotation_rate |
# or 
# | disp | volume | strain | stress | rotation |
#------------------------------------------------------------------------------
1   1   1   1   1
#------------------------------------------------------------------------------
#	OUTPUT FILES FOR AN ARBITRARY POINT DISLOCATION SOURCE
#               (for applications to earthquakes)
#	======================================================
# 1. selection (0 = not selected; 1 or 2 = selected), if (selection = 1), then
#    the 6 moment tensor elements [N*m]: Mxx, Myy, Mzz, Mxy, Myz, Mzx (x is
#    northward, y is eastward and z is downard); else if (selection = 2), then
#    Mis [N*m] = isotropic moment part = (MT+MN+MP)/3, Mcl = CLVD moment part
#    = (2/3)(MT+MP-2*MN), Mdc = double-couple moment part = MT-MN, Strike [deg],
#    Dip [deg] and Rake [deg].
#
#    Note: to use this option, the Green's functions above should be computed
#          (selection = 1) if they do not exist already.
# 2. switch for azimuth distribution of the stations (0 = uniform azimuth,
#    else = irregular azimuth angles)
# 3. list of the azimuth angles [deg] for all stations given above (if the
#    uniform azimuth is selected, then only one azimuth angle is required)
#
#------------------------------------------------------------------------------
#     Mis        Mcl        Mdc        Strike     Dip        Rake      File
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#     Mxx        Myy        Mzz        Mxy        Myz        Mzx       File
#------------------------------------------------------------------------------
0
0
0
#------------------------------------------------------------------------------
#
#	GLOBAL MODEL PARAMETERS (Note 5)
#	================================
# 1. switch for flat-earth-transform
# 2. gradient resolution [%] of vp, vs, and ro (density), if <= 0, then default
#    values (depending on wave length at cut-off frequency) will be used
#------------------------------------------------------------------------------
0
0  0  0
#------------------------------------------------------------------------------
#
#	                LAYERED EARTH MODEL
#       (SHALLOW SOURCE + UNIFORM DEEP SOURCE/RECEIVER STRUCTURE)
#	=========================================================
# 1. number of data lines of the layered model (source site) 
#------------------------------------------------------------------------------
25
#------------------------------------------------------------------------------
#
#	MULTILAYERED MODEL PARAMETERS (source site)
#	===========================================
# no  depth[km]    vp[km/s]    vs[km/s]   ro[g/cm^3]   qp      qs
#------------------------------------------------------------------------------
1  0.00  5.9000  3.4400  2.6700  927.34  599.99
2  15.53  5.9000  3.4400  2.6700  927.34  599.99
3  15.53  6.3000  3.6200  2.7400  927.34  599.99
4  25.26  6.3000  3.6200  2.7400  927.34  599.99
5  25.26  6.9000  3.8700  2.9100  927.34  599.99
6  35.00  6.9000  3.8700  2.9100  927.34  599.99
7  35.00  8.0574  4.4817  3.3509  927.34  599.99
8  49.50  8.0574  4.4817  3.3509  927.34  599.99
9  49.50  8.0545  4.4836  3.3864  927.34  599.99
10  51.50  8.0545  4.4836  3.3864  927.34  599.99
11  51.50  8.0539  4.4841  3.3939  927.34  599.99
12  53.00  8.0539  4.4841  3.3939  927.34  599.99
13  53.00  8.0493  4.4871  3.4492  927.34  599.99
14  77.50  8.0493  4.4871  3.4492  114.87  75.60
15  77.50  8.0477  4.4949  3.4642  114.87  75.60
16  120.00  8.0477  4.4949  3.4642  115.30  76.06
17  120.00  8.0816  4.5022  3.4128  115.30  76.06
18  142.50  8.0816  4.5022  3.4128  115.30  76.06
19  142.50  8.1438  4.5067  3.3849  115.30  76.06
20  165.00  8.1438  4.5067  3.3849  118.37  76.55
21  165.00  8.2064  4.5113  3.3593  118.37  76.55
22  187.50  8.2064  4.5113  3.3593  118.37  76.55
23  187.50  8.2692  4.5160  3.3359  118.37  76.55
24  210.00  8.2692  4.5160  3.3359  125.22  79.40
25  210.00  8.3309  4.5335  3.3313  125.22  79.40
#------------------------------------------------------------------------------
#
#	          LAYERED EARTH MODEL
#       (ONLY THE SHALLOW RECEIVER STRUCTURE)
#       =====================================
# 1. number of data lines of the layered model
#
#    Note: if the number = 0, then the receiver site is the same as the
#          source site, else different receiver-site structure is considered.
#          please be sure that the lowest interface of the receiver-site
#          structure given given below can be found within the source-site
#          structure, too.
#
#------------------------------------------------------------------------------
0
#------------------------------------------------------------------------------
#
#	MULTILAYERED MODEL PARAMETERS (shallow receiver-site structure)
#	===============================================================
# no  depth[km]    vp[km/s]    vs[km/s]   ro[g/cm^3]   qp      qs
#------------------------------------------------------------------------------
#---------------------------------end of all inputs----------------------------"""