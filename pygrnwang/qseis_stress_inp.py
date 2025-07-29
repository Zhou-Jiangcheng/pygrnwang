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
17.0
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
0.000
1  1
100
20.00 119.00
0.00 127.50 256
1 8.0
#------------------------------------------------------------------------------
#
#	WAVENUMBER INTEGRATION PARAMETERS
#	=================================
# 1. select slowness integration algorithm (0 = suggested for full wave-field
#    modelling; 1 or 2 = suggested when using a slowness window with narrow
#    taper range - a technique for suppressing space-domain aliasing);
# 2. 4 parameters for low and high slowness (Note 1) cut-offs [s/km] with
#    tapering: 0 < slw1 < slw2 defining cosine taper at the lower end, and 0 <
#    slw3 < slw4 defining the cosine taper at the higher end. default values
#    will be used in case of inconsistent input of the cut-offs (possibly with
#    much more computational effort);
# 3. parameter for sampling rate of the wavenumber integration (1 = sampled
#    with the spatial Nyquist frequency, 2 = sampled with twice higher than
#    the Nyquist, and so on: the larger this parameter, the smaller the space-
#    domain aliasing effect, but also the more computation effort); 
# 4. the factor for suppressing time domain aliasing (> 0 and <= 1) (Note 2).
#------------------------------------------------------------------------------
0
0.000  0.000  0.000  0.000
2.00
0.01
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
4   1
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
#    which will be appended by the program automatically: *.tz, *.tr, *.tt
#    and *.tv are for the vertical, radial, tangential, and volume change (for
#    hydrophones) components, respectively)
#------------------------------------------------------------------------------
#  explosion   strike-slip  dip-slip  clvd  single_f_v  single_f_h
#------------------------------------------------------------------------------
1    1    1    1    0    0
'ex' 'ss' 'ds' 'cl' 'fz' 'fh'
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
1
0  0  0
#------------------------------------------------------------------------------
#
#	                LAYERED EARTH MODEL
#       (SHALLOW SOURCE + UNIFORM DEEP SOURCE/RECEIVER STRUCTURE)
#	=========================================================
# 1. number of data lines of the layered model (source site) 
#------------------------------------------------------------------------------
43
#------------------------------------------------------------------------------
#
#	MULTILAYERED MODEL PARAMETERS (source site)
#	===========================================
# no  depth[km]    vp[km/s]    vs[km/s]   ro[g/cm^3]   qp      qs
#------------------------------------------------------------------------------
1  0.  5.8  3.46  2.6  1264.  600.
2  20.  5.8  3.46  2.6  1264.  600.
3  20.  6.5  3.85  2.9  1283.  600.
4  35.  6.5  3.85  2.9  1283.  600.
5  35.  8.04  4.48  3.58  1449.  600.
6  77.5  8.045  4.49  3.5  1445.  600.
7  77.5  8.045  4.49  3.5  180.6  75.
8  120.  8.05  4.5  3.427  180.  75.
9  120.  8.05  4.5  3.427  182.6  76.06
10  165.  8.175  4.509  3.371  188.7  76.55
11  210.  8.301  4.518  3.324  201.  79.4
12  210.  8.3  4.519  3.323  338.2  133.7
13  300.  8.628  4.679  3.401  353.6  138.7
14  410.  9.03  4.87  3.506  377.5  146.5
15  410.  9.36  5.08  3.929  414.1  162.7
16  660.  10.2  5.611  3.918  428.5  172.9
17  660.  10.79  5.965  4.24  1349.  549.5
18  764.  11.07  6.215  4.359  1276.  537.1
19  849.1  11.21  6.272  4.463  1263.  527.3
20  1038.  11.52  6.407  4.616  1230.  506.9
21  1227.  11.81  6.527  4.71  1198.  488.
22  1416.  12.08  6.636  4.812  1168.  470.4
23  1605.  12.33  6.735  4.909  1141.  454.4
24  1795.  12.55  6.827  5.004  1111.  437.7
25  1984.  12.78  6.913  5.096  1065.  415.6
26  2173.  13.  6.997  5.186  1037.  400.6
27  2362.  13.21  7.079  5.273  995.4  381.2
28  2551.  13.43  7.162  5.357  968.1  367.3
29  2740.  13.65  7.248  5.439  932.5  350.6
30  2740.  13.65  7.248  5.693  722.7  271.7
31  2790.  13.65  7.259  5.72  726.9  274.
32  2839.  13.66  7.27  5.746  725.1  274.
33  2892.  13.66  7.282  5.772  723.1  274.
34  2892.  7.972  0.  9.928  5.782E+04  0.
35  3344.  8.744  0.  10.59  5.782E+04  0.
36  3797.  9.338  0.  11.12  5.782E+04  0.
37  4249.  9.759  0.  11.55  5.782E+04  0.
38  4702.  10.09  0.  11.88  5.782E+04  0.
39  5154.  10.32  0.  12.14  5.782E+04  0.
40  5154.  11.04  3.505  12.7  632.8  85.03
41  5397.  11.12  3.567  12.82  619.2  85.03
42  5884.  11.23  3.647  12.97  605.  85.03
43  6371.  11.27  3.672  13.02  600.6  85.03
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
