# Prepare an Earth model

Use the same physical model for Green's functions, material normalization and
travel times. A model filename in one step does not override all other steps.

## Propagation input

Preprocessors accept `path_nd` containing six numeric columns:

```text
depth_km  vp_km_s  vs_km_s  density_g_cm3  Qp  Qs
```

Repeated depths represent discontinuities: rows above/below give the material
on each side. Preserve depth order and named boundaries such as `mantle`,
`outer-core` and `inner-core` where the TauP format uses them. Do not add
arbitrary prose, extra columns or headings: the conversion helpers distinguish
numeric rows from one-word boundary labels.

The examples use AK135 elastic velocities/density included in the package and
append illustrative constant `Qp=600` and `Qs=300`. This is **not** the
original AK135-F attenuation model. The QSEIS and EDGRN examples use the
first 24 numeric rows (ending at 809.5 km); spherical examples retain the
full Earth. No model download or developer-specific absolute path is needed.
Their shared helper is:

```{literalinclude} ../../examples/common.py
:language: python
:caption: Shared model, plotting and verification helpers
```

Resolve `path_nd` and `path_green` to absolute paths before launching a
backend. Some low-level calls change the current working directory, so
relative paths can fail after the first job starts.

## Extent and row selection

QSEIS uses a layered half-space with optional flat-Earth transformation.
Spherical SPGRN/QSSP require an appropriate radially layered Earth model.
A shallow crustal table is not automatically a full spherical model. Keep
the complete model in the spherical tutorials.

`earth_model_layer_num` controls the number of numeric model rows, not the
number of distinct geological layers. Leave it `None` for the supplied
model. If changed, inspect the generated input and confirm its declared
row count and final model boundary.

Zero shear velocity represents fluid. EDCMP moment normalization divides by
shear modulus and therefore requires a suitable solid source layer.
Strain-to-stress conversion also needs the intended receiver-layer material.

## Material values

`read_material_nd(model_name, depth)` returns
`[depth, vp, vs, density]` at or just below the requested depth. It recognizes
`ak135fc` or a **four-column** model filename, a different model-name
convention from TauP. Pass the generated `noQ.nd`, not the six-column
propagation input: `read_nd` defaults to four-column reshaping and cannot
infer the intended columns.

For solid material in these units,

$$
\mu=\rho v_s^2\,10^9,\qquad
\lambda=\rho(v_p^2-2v_s^2)\,10^9
$$

give Pa. Use source-layer material for EDCMP normalization and receiver-layer
material for local strain-to-stress conversion. At an interface, choose the
physically relevant side explicitly.

## Travel-time input

Preprocessing writes `noQ.nd` by removing attenuation columns. Java reads
this model directly; ObsPy builds an adjacent `.npz`, requiring a writable
model directory. `taup_create_npz_file` returns a path usable by the active
backend: `.nd` for Java, `.npz` for ObsPy.

Prepare a custom model before starting workers. Changing the propagation
model requires rebuilding travel-time tables as well; reuse flags do not
compare model contents. See [TauP](taup.md) and [resuming](parallel.md#resuming-a-calculation).

## Before scaling up

1. Run the original backend tutorial and retain its summary.
2. Substitute your model in a new output directory and inspect `grn.inp` or
   `spec.inp`.
3. Calculate a few P/S times at the intended depths and distances; check
   finite arrivals for the phases the reader will use.
4. Build one source/receiver pair and check units, shape and a known arrival
   or limiting case before expanding the grid.
