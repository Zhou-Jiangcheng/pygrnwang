# Scientific conventions

These conventions describe the Python readers in this repository. They differ
from some coordinate systems and component orders used by the Fortran programs.
Keep the model, source convention, source time function, units and time origin
with every exported array.

## Distances, models and units

| Quantity | Python convention |
| --- | --- |
| Source and receiver depth | km, positive downward |
| Epicentral distance | km; spherical conversion uses a radius of 6371 km |
| Azimuth `az_deg` | Degrees clockwise from north, from source to receiver |
| `sampling_interval`, `time_window` | Seconds |
| Reader `srate` | Samples per second (Hz) |
| Model `vp`, `vs`, density | km/s, km/s, g/cm³ in six-column `.nd` input |
| Model `Qp`, `Qs` | Dimensionless quality factors |
| Scalar moment and moment tensor | N m |
| Elastic moduli `lam`, `mu` | Pa |
| Source radius | km for SPGRN/QSSP; QSEIS2025 exposes dimensionless `source_radius_ratio` |
| Slowness | s/km for spherical backend input; SPGRN2020 native arrival tables use s/m; TauP ray parameter uses s/radian |

For moment-scaled dynamic synthetics, displacement is in m, velocity in m/s,
acceleration in m/s², strain is dimensionless, stress is in Pa, and rotation is
in radians. Rates add s⁻¹. A unit-moment calculation has the corresponding units
per N m. QSSP `gravitation` is a three-component acceleration; `gravimeter` is
scalar gravity change with **downward positive**, including ground acceleration
and the free-air-gradient effect.

QSEIS2025 `volume` is the fractional volume-change observable (volumetric
strain), not a volume in cubic metres. Its source-time-function handling follows
the strain family. The wrappers convert EDGRN/EDCMP inputs to SI internally;
do not pre-convert a `.nd` model to metres or kg/m³.

## Source mechanisms and amplitude

Dynamic readers call `check_convert_fm`:

| Input | Meaning |
| --- | --- |
| `[strike, dip, rake]` | Double couple, scalar moment 1 N m; angles in degrees |
| `[M0, strike, dip, rake]` | Double couple, scalar moment `M0` in N m |
| `[Mnn, Mne, Mnd, Mee, Med, Mdd]` | Six NED components in N m, used directly |
| `[M0, Mnn, Mne, Mnd, Mee, Med, Mdd]` | Tensor shape normalized and scaled to scalar moment `M0` |

The six-element order is **NN, NE, ND, EE, ED, DD**, not diagonal-first.
NED means north, east, down. The QSSP reader performs the conversion to
`[Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]` internally. Dynamic readers preserve the
amplitude of the converted tensor. Three angles alone produce a unit-moment
kernel, not a typical earthquake waveform amplitude.

Scalar moment is

$$
M_0 = \sqrt{\frac{M_{nn}^2+M_{ee}^2+M_{dd}^2
+2(M_{ne}^2+M_{nd}^2+M_{ed}^2)}{2}}.
$$

The seven-element form requires a nonzero tensor shape. An `M0` parameter
does not interpret its input as moment magnitude.

### EDCMP normalization

`seek_edcmp2` and its bulk counterpart normalize the mechanism to unit scalar
moment, even if an amplitude was supplied. `check_convert_pure_dp=True` also
projects the shape onto a double couple. With `False`, the normalized tensor
is synthesized from the five available deviatoric basis sources; this is not
a full isotropic-source library.

The basis dislocations use unit slip and unit area. By default,
`times_mu=False` divides the result by source-layer shear modulus, giving a
kernel per N m. Multiply this default output by `M0` to obtain physical static
deformation. Supplying `[M0, strike, dip, rake]` alone does **not** scale EDCMP
output.

`times_mu=True` leaves the raw unit-slip, unit-area normalization.
`area_km_sq` always multiplies by `area_km_sq * 1e6`, whichever normalization
was selected. Thus `times_mu=True, area_km_sq=A` yields the unit-slip response
for area `A`; multiply by slip in metres. Do not also multiply that result by
seismic moment.

The modulus lookup uses `rho * vs**2 * 1e9` for model units. Pass the same
material model used to build the library. The lookup accepts `ak135fc` or a
four-column model filename such as `noQ.nd`. Do not pass the six-column
propagation input to the material reader. `ak135` is a TauP model name,
not a built-in material-model option.

## Vector components

Main dynamic `seek_*` readers return `(n_components, n_samples)`, including
`(1, n_samples)` for a scalar observable. For vectors, `rotate=True` returns
**E, N, Z**, with Z upward; it does not return NED. `rotate=False` returns
**R, T, Z**. R points outward from the source and T points to the left of R
viewed from above. At azimuth 0°, R is north and T is west.

The implemented transform is

$$
E=R\sin a-T\cos a,\qquad N=R\cos a+T\sin a,\qquad Z=Z.
$$

This order also applies to QSEIS2025/QSSP rotation observables. EDCMP
displacement is a one-dimensional vector of length 3; its tilt has length 2,
ordered E/N when rotated and R/T otherwise.

## Six-component tensors

Rotated tensors use **EE, EN, EZ, NN, NZ, ZZ**. Off-diagonal strains are tensor
strains, not engineering shear strains; stress conversion uses
`2 * mu * strain_en`.

| Reader | `rotate=True` | `rotate=False` |
| --- | --- | --- |
| QSEIS2025 strain/stress and rates | EE, EN, EZ, NN, NZ, ZZ | North-reference ENZ basis |
| QSSP2020 strain/stress and rates | EE, EN, EZ, NN, NZ, ZZ | North-reference ENZ basis |
| EDCMP strain/stress | EE, EN, EZ, NN, NZ, ZZ | RR, RT, RZ, TT, TZ, ZZ |
| QSEIS06 main reader; SPGRN2012/2020 | No tensor observable | No tensor observable |

The **north-reference basis** is the synthesized tensor before the final
azimuthal rotation. It differs from the vector reader's RTZ output.
Expressed in the reader's local RTZ axes, its entries correspond to
`[TT, -TR, -TZ, RR, RZ, ZZ]`. In QSEIS2025 source-file notation the actual
assembled array is `[s_tt, s_rt, -s_zt, s_rr, -s_zr, s_zz]`, with `e_*`
instead for strain. Fortran uses downward z and north-to-east transverse t.

Use `rotate=True` for tensor comparisons between backends. Do not apply a
vector rotation to six components or assume `rotate=False` has one universal
tensor order.

The separate QSEIS06 finite-difference readers construct strain-rate and
stress-rate from extra receiver-depth and distance samples. They have no
`rotate` argument and use their own azimuthal transformation. Their accuracy
and signs need a dedicated derivative/convergence check; use the direct
QSEIS2025 tensor workflow for the introductory tensor example.

## Source time functions

| Backend | Parameter and unit | Interpretation |
| --- | --- | --- |
| QSEIS06/2025 | `wavelet_duration`, **number of samples** | Integer duration; multiply by `sampling_interval` to interpret in seconds |
| SPGRN2012/2020 | `source_duration`, seconds | Squared half-sinusoid duration |
| QSSP2020 | `source_duration`, seconds | Squared half-sinusoid moment-rate duration |
| EDGRN/EDCMP | None | Static response |

QSEIS `wavelet_type=1` selects a normalized squared half-sinusoid approximating
a delta impulse; stored vector kernels represent velocity. Type 2 selects its
integral, a tapered Heaviside; stored vector kernels represent displacement.
The reader integrates/differentiates to obtain the requested observable.
The same distinction applies to rate/non-rate strain, stress, volume and
rotation kernels in QSEIS2025. A nonpositive duration requests the Fortran
default of two samples, not an infinitely short physical source. The writer
formats duration as an integer.

SPGRN stores velocity kernels and its reader integrates displacement or
differentiates acceleration. QSSP writes each selected observable directly.
For the spherical examples, a duration `T=64` s describes a normalized
squared half-sinusoid on 0 to `T`, with its centroid at 32 s. The examples
do not apply a separate centroid shift.

There is an implementation difference despite the shared physical source
definition: SPGRN2012's old wavelet routine omits the imaginary frequency
used in numerical damping, whereas SPGRN2020 and QSSP2020 include it.
For the examples' 4096 s FFT period and 0.01 anti-aliasing factor, the
effective SPGRN2012 source-pulse area after damping correction is about
3.67% larger. This value depends on the duration and damping settings;
it is not a universal amplitude conversion. See the
[backend comparison](guides/backend-comparison.md) before interpreting
small inter-backend amplitude differences.

The QSEIS regional examples use `wavelet_type=0` with 1024 custom
moment-rate nodes spanning 64 s. For target rate
`r(t) = (2/64) * sin(pi*t/64)**2` on 0–64 s, the written samples are
`r(t) * exp(2*pi*fi*t)`, where
`fi = log(0.01) / (2*pi*4092)`. This compensates for the solver's
real-frequency wavelet transform and subsequent damping correction.
The input samples have area about 0.9647094 and must not be renormalized:
the effective physical pulse has unit area and centroid 32 s.
The verified time-domain and spectral relative L2 errors are below
`1e-5` against the target. See the
[STF construction and checks](guides/backend-comparison.md#matching-the-effective-source-time-function).

Custom QSEIS wavelets require a sample block in the low-level input;
the high-level preprocessor has no custom-array parameter. The regional
example helper installs this block after preprocessing. The ordinary
readers do not infer the normalization or rate/non-rate meaning of type 0.
For this moment-rate pulse, the examples explicitly read `velo`,
`strain_rate` or `stress_rate` and integrate once with `cumsum * dt`.
Calling the generic reader with `output_type="disp"` on this type-0
library does not automatically integrate velocity. The default
near-distance tutorials retain their built-in type-2 pulse.

## Time origin, reduction and arrivals

At the native sampling interval `dt`, unadjusted sample `i` represents
`t_start + i * dt`:

| Backend | Native `t_start` relative to source origin |
| --- | --- |
| QSEIS06/2025 | `grn_dist / time_reduction_velo`, or zero for zero reduction velocity |
| SPGRN2012 | Nearest integer second to `t0 + distance / v0` |
| SPGRN2020 | Nearest integer second to direct-P onset minus `green_before_p` |
| QSSP2020 | `time_reduction`, in seconds |

SPGRN native binary record headers contain the stored start time.
Use those values for a strict source-origin comparison: metadata formulas
or fractional P onsets alone can miss integer-second rounding.
SPGRN2012 reader alignment uses requested distance and `v0`; use exact grid
points for simple timing comparisons. Negative QSSP reduction means before
source origin, not before P. QSEIS reduction velocity is in km/s.

`before_p` selects seconds preceding the library P onset in the returned
array. `pad_zeros=True` shifts a reduced-time trace toward source-origin time
and fills exposed samples with zeros. Do not combine these options. They
preserve array length and can discard samples at the opposite end.

`shift=True` recomputes P/S times at the requested location and applies
piecewise interpolation between arrival regions. It approximates a waveform;
it does not calculate a new Green's function. Validate against a direct
calculation before phase-sensitive work or crossing phase branches. Finite
P and S times are required.

With `only_seismograms=False`, main dynamic readers return
`(seismograms, tpts_table, first_p, first_s, source_depth, receiver_depth, distance)`.
The last three values identify nearest library nodes even with waveform
interpolation. `first_p` and `first_s` remain `None` when `shift=False`,
including with `before_p`. QSEIS/QSSP also leave `tpts_table=None` unless a
timing operation needs it. SPGRN reads its table on every query. A missing
physical arrival in TauP is `NaN`, distinct from `None`.

See [reading](guides/reading.md) before assigning plot time axes. These
conventions were checked against current Python assembly/rotation code and
Fortran input/output routines. Inconsistencies are recorded in
[known limitations](guides/troubleshooting.md#known-implementation-limitations).

## Convergence and reproducibility

Choose the shortest period and latest arrival your study needs, then
converge frequency/slowness or harmonic cutoffs, time window, wavenumber
sampling, source duration and model discretization. A small tutorial
demonstrates the workflow; it does not establish numerical accuracy for
another source, distance or frequency band. Save the exact model, generated
inputs, package version, grid and processing settings with published results.

For spherical backends, equal positive slowness limits do not guarantee
equal low-frequency harmonic content. In SPGRN2020, `max_slowness=0`
selects a separate full-wavefield branch with automatic slowness selection.
In QSSP2020, `min_harmonic` controls the low-frequency baseline of an
upper cutoff, while `max_harmonic` also affects spatial differential
transformation. Neither parameter describes a band that drops all lower
degrees. The [controlled example comparison](guides/backend-comparison.md)
shows why file, shape and finite-value checks alone cannot establish
waveform accuracy.
