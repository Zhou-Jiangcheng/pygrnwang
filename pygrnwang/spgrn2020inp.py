s = """# This is the input file of FORTRAN77 program "spgrn2020" for calculating
# synthetic Green's functions of a self-gravitating, spherically symmetric,
# isotropic and viscoelastic earth.
#
# by
# Rongjiang  Wang <wang@gfz-potsdam.de>
# Helmholtz-Centre Potsdam
# GFZ German Reseach Centre for Geosciences
# Telegrafenberg, D-14473 Potsdam, Germany
#
# Last modified: Potsdam, April, 2020
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# If not specified, SI Unit System is used overall!
#
# Coordinate systems:
# spherical (r,t,p) with r = radial,
#                        t = co-latitude,
#                        p = east longitude.
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#	UNIFORM RECEIVER DEPTH
#	======================
# 1. uniform receiver depth [km]
#-------------------------------------------------------------------------------------------
0.00
#-------------------------------------------------------------------------------------------
#
#	TIME (FREQUENCY) SAMPLING
#	=========================
# 1. time window [s] and sampling interval [s]
# 2. max. frequency [Hz] of Green's functions
# 3. max. slowness [s/km] of Green's functions (if <= 0.0, then complete wave field
#    including static deformation will be calculated)
# 4. anti-aliasing factor (> 0 & < 1), if it is <= 0 or >= 1/e (~ 0.4), then
#    default value of 1/e is used (e.g., 0.1 = alias phases will be suppressed
#    to 10% of their original amplitude)
#
#    Note: The computation effort increases linearly the time window and
#          quadratically with the cut-off frequency.
#-------------------------------------------------------------------------------------------
4095.75  0.25
1.00
0.00
0.1
#-------------------------------------------------------------------------------------------
#
#	SELF-GRAVITATING EFFECT
#	=======================
# 1. the critical frequency [Hz] and the critical harmonic degree, below which
#    the self-gravitating effect should be included # 0.03 300
#-------------------------------------------------------------------------------------------
0.0  0
#-------------------------------------------------------------------------------------------
#
#	WAVE TYPES
#	==========
# 1. selection (1/0 = yes/no) of speroidal modes (P-SV waves), selection of toroidal modes
#    (SH waves)
#-------------------------------------------------------------------------------------------
1     1
#-------------------------------------------------------------------------------------------
#	SPECTRAL GREEN'S FUNCTIONS
#	==========================
# 1. directory for spectral Green's functions
# 2. number of discrete source depths
# 3. list of source depths [km], source radius [km], the respective file names of Green's
#    functions and switch number (0/1) (0 = spectra of Green's functions for this depth
#    exist already, 1 = need to be calculated or updated.
#    Note: update is required if any of the above parameters is changed)
#-------------------------------------------------------------------------------------------
''
1
10.0  0.00  'grn_d10.0'  1
#-------------------------------------------------------------------------------------------
#	SPACE-TIME DOMAIN GREEN'S FUNCTIONS
#	===================================
# 1. directory for space-time domain Green's functions
# 2. header file of Green function database
# 3  file for tabulating P wave arrivals
# 4  file for tabulating S wave arrivals
# 5. output time window [sec] and sampling interval [s]
# 6. start time of seismograms in reference to the direct P wave onset (< 0): t0[s]
# 7. source wavelet (square half-sinusoid) duration [s]
# 8. selection of order of Butterworth band-pass filter (if <= 0, then no filtering), lower
#    and upper cutoff frequencies [hz]
# 9. epicentral distance sampling used for Green's functions: d1, d2, delta_1, delta_2 [km],
#    where delta_1 and delta_2 are the smallest and largest intervals near d1 and d2,
#    respectively (sampling interval linearly increasing with distance)
#-------------------------------------------------------------------------------------------
''
'GreenInfo10.0.dat'
'tptable.dat'
'tstable.dat'
4095.5  0.5
-120.0
0.00
0  0.005  0.500
3300  11150  10.0  10.0
#--------------------------------------------------------------------------------------------------------
#
#   MULTI-LAYERED EARTH MODEL (AK135)
#   =================================
# 1. number of data lines of the layered model and selection for including
#    the physical dispersion according to Kamamori & Anderson (1977)
#--------------------------------------------------------------------------------------------------------
43  0
#--------------------------------------------------------------------------------------------------------
#
#   MODEL PARAMETERS
#   ================
# no    depth[km] vp[km/s]  vs[km/s]   ro[g/cm^3] qp         qs
#--------------------------------------------------------------------------------------------------------
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
#---------------------------------end of all inputs------------------------------------------------------"""
