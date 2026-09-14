# EDGRN → EDCMP: static deformation

EDGRN computes static Green's tables in a layered elastic half-space.
EDCMP combines those tables for dislocation sources. The Python library
workflow builds five basis dislocations and synthesizes static
displacement, strain, stress or tilt at a query point.

## Complete layered calculation

```console
python examples/edgrn_edcmp.py
```

The script runs four necessary steps: EDGRN preprocessing/calculation,
EDCMP preprocessing/calculation using the same tables, format conversion,
and displacement reading with explicit seismic-moment scaling.

```{literalinclude} ../../examples/edgrn_edcmp.py
:language: python
:caption: Complete layered EDGRN and EDCMP workflow
```

```{figure} ../_static/examples/edgrn_edcmp.png
:alt: Static east, north and up displacement versus epicentral distance.

Static displacement in metres. There is no time axis or sampling interval.
```

The model uses the first 24 numeric rows. EDGRN requires at least two
source depths, so the library contains 10 and 11 km and the example
queries 10 km. Both stages use distances 0, 30, 60, 90 and 120 km, while
plots query the interior 30/60/90 km points. The extra margin avoids
finite-dislocation geometry crossing the table boundary.

The default output is `examples/output/edgrn_edcmp/` with `disp.npz`,
`disp.png` and `summary.json`. The saved array has shape `(3, 3)`:
one E/N/up vector per queried distance.

## Parameters and two-stage consistency

Keep source depth range/increment, receiver depth list, distance
range/increment and root directory consistent between EDGRN and EDCMP.
The wrappers accept km and convert the backend inputs to metres.
`wavenumber_sampling_rate` controls EDGRN numerical integration;
converge it for your geometry and material model.

`layered=True` makes EDCMP use EDGRN's `edgrn.ss`, `edgrn.ds` and
`edgrn.cl`. The low-level EDCMP writer also has a homogeneous half-space
mode with `layered=False` and `lam`/`mu` in Pa. The complete tutorial
uses the layered mode and the model-based normalization.

The four `output_observables` positions are displacement, strain,
stress and tilt. The example selects `(1, 0, 0, 0)`. To build another
quantity, enable its flag before calculation and request its corresponding
`output_type` from the reader.

## Normalize the static response

EDCMP readers normalize a supplied mechanism to scalar moment 1, regardless
of whether the mechanism contains an explicit `M0`. The example requests
the default `times_mu=False` kernel and then multiplies by `10^15 N m`.

The material lookup must receive the generated **four-column `noQ.nd`**,
not the six-column propagation model. Passing the latter can produce an
incorrect modulus through the low-level reshape without a useful error.
The script supplies the correct material file.

Alternatively, `times_mu=True` retains unit-slip, unit-area basis values.
Combined with `area_km_sq=A` and an external slip factor in metres, it
gives the corresponding dislocation response. The bulk reader also accepts
`slip_m_arr`. Do not multiply that area/slip-scaled result by `M0` again.
See the full [normalization contract](../conventions.md#edcmp-normalization).

`check_convert_pure_dp=True` projects source shape to a double couple.
Disabling it preserves a normalized tensor shape within the five-basis
representation; it does not restore input magnitude or create an
isotropic-source basis.

## Files and result interpretation

EDGRN writes under `edgrn2/<receiver>/`; EDCMP writes five source-basis
directories under `edcmp2/<source>/<receiver>/`. The sequential EDCMP
builder does not call the bulk converter automatically, so the script
explicitly invokes `convert_pd2bin_edcmp2_all`. The converted whole-grid
files are also needed for `seek_edcmp2_bulk`.

Single queries return vectors: length 3 for displacement, 6 for
strain/stress and 2 for tilt. Rotated strain/stress order is
EE/EN/EU/NN/NU/UU; unrotated order is RR/RT/RZ/TT/TZ/ZZ.
Tensor strain is dimensionless, stress is Pa and tilt is rad after the
appropriate amplitude scaling.

These are static results, not the final sample of a dynamic waveform
with an arbitrary time window. A comparison with dynamic late-time
response requires matching source normalization, half-space physics,
converged low-frequency content and a sufficiently long trace.
