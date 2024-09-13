s = """# This is the input file of FORTRAN77 program "qssp2020" for calculating
# synthetic seismograms of a self-gravitating, spherically symmetric,
# isotropic and viscoelastic earth.
#
# by
# Rongjiang Wang <wang@gfz-potsdam.de>
# Helmholtz-Centre Potsdam
# GFZ German Reseach Centre for Geosciences
# Telegrafenberg, D-14473 Potsdam, Germany
#
# Last modified: Potsdam, April 2020
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# If not specified, SI Unit System is used overall!
#
# Coordinate systems:
# spherical (r,t,p) with r = radial, t = co-latitude, and p = east longitude.
# local cartesian (e,n,z) with e = east, n = north, and z = vertical (upwards positve).
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#   UNIFORM RECEIVER DEPTH
#   ======================
# 1. uniform receiver depth [km]
#--------------------------------------------------------------------------------------------------------
    0.00
#--------------------------------------------------------------------------------------------------------
#
#   SPACE-TIME SAMPLING PARAMETERS
#   ==============================
# 1. time window [sec], sampling interval [sec]
# 2. max. frequency [Hz] of Green's functions
# 3. max. slowness [s/km] of Green's functions
#    Note: if the near-field static displacement is desired, the maximum slowness should not
#          be smaller than the S wave slowness in the receiver layer
# 4. anti-aliasing factor (> 0 & < 1), if it is <= 0 or >= 1/e (~ 0.4), then default value of
#    1/e is used (e.g., 0.1 = alias phases will be suppressed to 10% of their original
#    amplitude)
# 5. switch (1/0 = yes/no) of turning-point filter, the range (d1, d2) of max. penetration
#    depth [km] (d1 is meaningless if it is smaller than the receiver/source depth, and
#    d2 is meaningless if it is equal to or larger than the earth radius)
#
#    Note: The turning-point filter (Line 5) works only for the extended QSSP code (e.g.,
#          qssp2016). if this filter is selected, all phases with the turning point
#          shallower than d1 or deeper than d2 will be filtered.
#
# 6. Earth radius [km], switch of free-surface-reflection filter (1/0 = with/without free
#    surface reflection)
#
#    Note: The free-surface-reflection filter (Line 6) works only for the extended QSSP
#          code (e.g., qssp2016). if this filter is selected, all phases with the turning
#          point shallower than d1 or deeper than d2 will be filtered.
#--------------------------------------------------------------------------------------------------------
    255.0   1.00
    0.05
    0.40
    0.01
	0   2891.5   6371.0
	6371.0    1
#--------------------------------------------------------------------------------------------------------
#
#   SELF-GRAVITATING EFFECT
#   =======================
# 1. the critical frequency [Hz] and the critical harmonic degree, below which
#    the self-gravitating effect should be included
#--------------------------------------------------------------------------------------------------------
    0.05     500
#--------------------------------------------------------------------------------------------------------
#
#   WAVE TYPES
#   ==========
# 1. selection (1/0 = yes/no) of spheroidal modes (P-SV waves), selection of toroidal modes
#    (SH waves), minimum and maximum cutoff harmonic degrees
#    Note: if the near-field static displacement is desired, the minimum cutoff harmonic
#          degree should not be smaller than, e.g., 2000.
#--------------------------------------------------------------------------------------------------------
    1     1    20000  20000
#--------------------------------------------------------------------------------------------------------
#   GREEN'S FUNCTION FILES
#   ======================
# 1. number of discrete source depths, estimated radius of each source patch [km] and
#    directory for Green's functions
# 2. list of the source depths [km], the respective file names of the Green's functions
#    (spectra) and the switch number (0/1) (0 = do not calculate this Green's function because
#    it exists already, 1 = calculate or update this Green's function).
#    Note: Green's functions need to be recalculated if any of the above parameters is changed.
#--------------------------------------------------------------------------------------------------------
    1     0.0  '/home/zjc/Desktop/test_qssp_gf/'
    0.00   'Green_0km'    0
#--------------------------------------------------------------------------------------------------------
#
#   MULTI-EVENT SOURCE PARAMETERS
#   =============================
# 1. number of discrete point sources and selection of the source data format
#    (1, 2 or 3)
# 2. list of the multi-event sources
#
#    Format 1 (full moment tensor):
#    Unit     Mrr  Mtt  Mpp  Mrt  Mrp  Mtp  Lat   Lon   Depth  T_origin T_rise
#    [Nm]                                   [deg] [deg] [km]   [sec]    [sec]
#
#    Format 2 (double couple):
#    Unit   Strike    Dip       Rake      Lat   Lon   Depth  T_origin T_rise
#    [Nm]   [deg]     [deg]     [deg]     [deg] [deg] [km]   [sec]    [sec]
#
#    Format 3 (single force):
#    Unit      Feast    Fnorth  Fvertical    Lat   Lon   Depth  T_origin T_rise
#    [N]                                     [deg] [deg] [km]   [sec]    [sec]
#
#    Note: for each point source, the default moment (force) rate time function is used, defined by a
#          squared half-period (T_rise) sinusoid starting at T_origin.
#--------------------------------------------------------------------------------------------------------
    1      1
    1  1.0 0.0 0.0 0.0 0.0 0.0  0.0 0.0 10.0   0.0 0.0
#--------------------------------------------------------------------------------------------------------
#
#   RECEIVER PARAMETERS
#   ===================
# 1. select output observables (1/0 = yes/no)
#    Note: "gravitation" here means (space-based) garvitational acceleration vector, while "gravity"
#          means the gravity change (downwards positive) measured by a gravimeter. The latter includes
#          the effect due to free-air gradient and ground acceleration.
# 2. output file name
# 3. output time window [sec] (<= Green's function time window)
# 4. selection of order of Butterworth bandpass filter (if <= 0, then no filtering), lower
#    and upper corner frequencies (<= cut-off frequency defined above)
# 5. lower and upper slowness cut-off [s/km] (slowness band-pass filter)
# 6. number of receiver
# 7. list of the station parameters
#    Format:
#    Lat     Lon    Name     Time_reduction
#    [deg]   [deg]           [sec]
#    (Note: Time_reduction = start time of the time window)
#--------------------------------------------------------------------------------------------------------------------
# disp | velo | acce | strain | strain_rate | stress | stress_rate | rotation | rotation_rate |gravitation |  gravity
#--------------------------------------------------------------------------------------------------------------------
   1      1      1      1              1        1        1            0        0            0      0
  '/home/zjc/Desktop/test_qssp_gf_output/test'
  256.0
  0   0.05   0.25
  0.00   0.40
  1
#--------------------------------------------------------------------------------------------------------
#
#   MULTI-LAYERED EARTH MODEL
#   =========================
# 1. number of data lines of the layered model and selection for including
#    the physical dispersion according to Kamamori & Anderson (1977)
#--------------------------------------------------------------------------------------------------------
    10     0
#--------------------------------------------------------------------------------------------------------
#
#   MODEL PARAMETERS
#   ================
# no   depth[km]   vp[km/s]   vs[km/s]   ro[g/cm^3] qp         qs
#--------------------------------------------------------------------------------------------------------
  1      0.000     5.8000     3.2000     2.6000    1478.30     599.99
  2      3.300     5.8000     3.2000     2.6000    1478.30     599.99
  3      3.300     5.8000     3.2000     2.6000    1478.30     599.99
  4     10.000     5.8000     3.2000     2.6000    1478.30     599.99
  5     10.000     6.8000     3.9000     2.9200    1368.02     599.99
  6     18.000     6.8000     3.9000     2.9200    1368.02     599.99
  7     18.000     8.0355     4.4839     3.6410     950.50     394.62
  8     43.000     8.0379     4.4856     3.5801     972.77     403.93
  9     80.000     8.0400     4.4800     3.5020    1008.71     417.59
 10     80.000     8.0450     4.4900     3.5020     182.03      75.60
#---------------------------------end of all inputs------------------------------------------------------"""
