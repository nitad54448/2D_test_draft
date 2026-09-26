# Thermoelectric Lab — 2D browser edition

A standalone HTML/CSS/JavaScript workbench for nonlinear planar thermoelectric
transport. The model solves spatial temperature, voltage and current fields in
multiple materials, including Seebeck, Joule, Peltier and Thomson effects. Periodic
runs return DC and the 1ω, 2ω and 3ω components.

## Start

Extract the ZIP and open **index.html**. No runtime packages or internet connection
are needed. All JavaScript, styles and worker source are embedded in that file.

To use Python's local web server, open a terminal in this extracted folder:

```bash
python -m http.server 8001 --bind 127.0.0.1
```

Open **http://127.0.0.1:8001**. Use `python3` if that is your Python command.
Port 8001 avoids interfering with a 1D app already running on port 8000. Stop the
server with Ctrl+C. Python only serves static files; the solver runs in JavaScript.

Keep this version in its own folder. It does not overwrite the 1D application.
All files in the ZIP are complete.

## The five tabs

| Tab | Controls / output |
|---|---|
| Geometry | Width, height, out-of-plane depth, x/y cell counts, paintable material grid |
| Materials | Material reference properties and temperature-dependent resistivity/Seebeck laws |
| Boundaries | Two equipotential electrical contacts, their positions/coverage, four thermal sides |
| Solver | Steady or periodic mode, drive frequency, steps per period, maximum cycles |
| Results | 2D heatmaps, harmonic amplitude/phase/real part, current arrows, movable temperature probe, terminal spectrum and exports |

Tabs retain their values. Arrow keys, Home and End navigate the tab list.

Suggested first runs:

1. **Open circuit · Seebeck DC** — checks the sign and value of Seebeck voltage.
2. **Narrow contact · current spreading** — a genuinely 2D stationary example.
3. **Homogeneous · Joule 2ω** — time-periodic Joule heating.
4. **Nonlinear resistance · 3ω** — voltage mixing through temperature-dependent resistivity.
5. **Cu / BiTe · layered** — coupled multilayer response.

Preset properties are illustrative, not certified material data.

## Geometry and materials

The domain is a filled rectangle in x/y, extruded through a constant physical
depth. x increases to the right; y increases upward. GUI dimensions are in mm;
the model configuration uses meters. Depth converts current density into amperes
and heat flux into watts. This is planar 2D, not axisymmetric geometry.

Click a palette material, then click/drag on the grid to paint its cells. Add
materials in the Materials tab. Painting is disabled while the worker is running.
The grid can represent arbitrary *grid-aligned* material regions, not holes,
unmeshed voids or curved boundaries. Each cell must have a conducting material.

After changing dimensions or cell counts, click **Apply dimensions / remesh**.
Remeshing resamples the old cell map by cell-center location. Inspect small regions
after remeshing; they can disappear on a coarse grid. All interfaces follow grid
lines. Electrode endpoints are resolved on boundary nodes and effectively snap
inward to nodes within the selected range.

The editor uses:

```
rhoe(T) = [1 + beta*(T - 300)] / sigma300
alpha(T) = alpha300 + alphaSlope*(T - 300)
```

Density, heat capacity and thermal conductivity are constant in the GUI. Seebeck
is entered in µV/K and its slope in µV/K². Temperatures are absolute K. The solver
requires positive finite density, heat capacity, electrical conductivity and
thermal conductivity. Alpha may have either sign.

The JS material class also accepts functions of T or serialized `linear` and
`inverseLinear` descriptors. Custom functions must be defined in source; they
cannot be transferred through worker messages. The GUI imports scalar reference
properties and its built-in linear laws; it does not evaluate text as code.

## Electrical contacts

Select a source and a sink side: left, right, bottom or top. Each contact is an
equipotential segment occupying a percentage range on that side. Ranges go:

- bottom to top on left/right sides;
- left to right on bottom/top sides.

White lines show the source, pink lines the sink. Each contact requires at least
two boundary nodes; overlapping contacts, including a shared corner, are rejected.
All other exterior nodes have zero normal electrical current.

The source potential is zero. Control modes are:

- **Voltage:** prescribe V(sink). V(source)=0.
- **Total current:** solve for the sink potential so the total current entering
  through the source equals the requested amperes. Local current density is not
  prescribed and can vary strongly over the contact.
- **Open circuit:** net current through each contact is zero, but each contact
  remains equipotential. Internal thermoelectric circulating currents can still
  exist in 2D. This does not impose J=0 in every material.

Positive total current enters at the source and exits at the sink. Reported
terminal voltage is **V(sink)−V(source)**. A passive isothermal resistor therefore
has negative terminal voltage under positive current.

Electrical excitation is `bias + peak*cos(2*pi*f*t + phase)`; phase is in degrees.
Steady mode uses only the DC bias. Periodic mode retains both bias and AC.
There is no external circuit capacitance or inductance.

## Thermal boundaries

Each of the four sides can use:

| Condition | Meaning |
|---|---|
| Temperature | Prescribed absolute temperature in K |
| Flux | Prescribed **total outward** heat flux in W/m² |
| Convection | Total outward heat flux h(T−Tambient), h in W/m² K |

Zero total flux is adiabatic, including Peltier transport. A fixed-temperature
boundary absorbs/provides the required reaction heat. Adjacent fixed-temperature
sides must agree at their common corner; contradictory values are rejected.
Temperature and flux boundary values may contain DC and AC components. The GUI
uses phase zero for thermal signals; the JS boundary descriptors can also specify
phase. Thermal electrodes and electrical contacts are independent settings.
There are no losses through the front/back faces associated with the extrusion
depth. Include those only after explicitly extending the model.

## Equations and conservative spatial model

```
J = -sigma * (grad(V) + alpha*grad(T))
div(J) = 0
q = alpha*T*J - k*grad(T)
rho*Cp*dT/dt = -div(q) - J.dot(grad(V))
```

The model uses scalar isotropic properties. Inside a smooth material, expansion
yields Joule heating J²/sigma and Thomson heating −T(dalpha/dT)J·grad(T).
Material jumps produce Peltier transport through the total heat flux. Those terms
must not be added again as separate interface or volumetric sources.

The implementation is a nodal conservative control-volume network on a rectangular
grid. Each material cell contributes four links between adjacent corner nodes:
two x links and two y links, each with half the cell's transverse face area.
Neighboring material cells can contribute separate parallel links along the same
node pair. This retains each side's material rather than averaging unrelated
materials into one nodal coefficient. Temperature and potential are continuous
at shared interface nodes; contacts between materials are ideal.

For a link a→b in one material:

```
g = Ahalf / (length * mean(rhoe(Ta), rhoe(Tb)))
Iab = g * [Va - Vb - alphaMean*(Tb - Ta)]
Qab = alphaMean*(Ta+Tb)/2 * Iab - kMean*Ahalf/length*(Tb-Ta)
Pab = Iab*(Va-Vb)
```

The thermal node balance receives −Qab at a, +Qab at b and half Pab at each end.
Each cell contributes one quarter of its local heat capacity to each corner.
Internal link fluxes cancel in global balances. An interface node receives
capacity contributions from all its neighboring materials.

These are finite-grid approximations. Small inclusions, sharp corners, high
property contrasts and contact singularities need spatial refinement. The app
is not a mesh-independent description of point contacts.

## Electrical and nonlinear solution

For fixed T, the electric graph problem is linear. Voltage mode solves it once.
Current/open-circuit mode combines a Seebeck-driven zero-terminal-voltage solution
with a unit-terminal-voltage solution, then chooses the terminal voltage satisfying
the total-current constraint. Current is a spatial field, not the uniform scalar
used by the 1D reduction.

Electric and thermal linear systems are symmetric graph operators solved using
Jacobi-preconditioned, matrix-free conjugate gradients. Fixed boundary nodes are
eliminated. There is no dense numerical Jacobian.

Nonlinear thermal/electric coupling uses damped Picard iterations. Each iteration
re-evaluates material laws, solves voltage/current, assembles thermal transport
and solves temperature. Default relaxation is 0.85, update tolerance 2e-9 K and
maximum 100 nonlinear iterations. Failure raises an explicit error. Picard can
fail on difficult strongly nonlinear cases; this version does not silently
switch algorithms or accept an unfinished iterate.

Periodic runs use fixed-step BDF2 with one backward-Euler startup step. Consecutive
full temperature cycles are compared using absolute tolerance 2e-7 K plus relative
tolerance 1e-10 of the absolute temperature. At least three cycles are required.
Default GUI settings: 128 steps per period and a 100-cycle budget. All serialized
boundary signals share the specified drive frequency, so no arbitrary nonperiodic
function is accepted as a boundary signal.

**A converged cycle does not guarantee accurate small harmonics.** Repeat with
twice as many cells and twice as many steps. For weak 3ω responses, start with
256 or 512 steps and compare refinements. Tighten nonlinear/cycle tolerances through
the API when necessary. Strong heating may have no stable periodic state.

The app caps the mesh at 1600 nodes. Start around 12×8 cells. Time and memory costs
increase with cell count, time samples, nonlinear iterations and cycles. The worker
keeps the page responsive; Stop terminates the calculation. No result is marked
converged after a solver failure.

## Harmonics and field display

```
u(t) = U0 + Re[sum(Un * exp(i*n*omega*t))]
```

All phasors use **peak**, not RMS, amplitudes. `{re, im}` stores each complex
coefficient. A positive sine has a negative imaginary phasor.

- Temperature and voltage are nodal. Their heatmaps display an average of the
  four corner values per cell; phase uses the average complex phasor before
  taking its angle. The color scale spans the underlying nodal values.
- Jx, Jy and heat-flux components are reconstructed at cell centers from the two
  corresponding half-face contributions. Cell-average reconstructions are for
  plotting; current conservation is enforced on the original link network.
- The current-magnitude map uses `sqrt(|Jx_n|² + |Jy_n|²)` for harmonic n. It is a
  vector phasor norm, not a Fourier coefficient of the instantaneous scalar |J|.
  Its phase is not defined, so the phase control is disabled for this map.
- DC displays signed real values, except |J|. At higher harmonics select amplitude,
  phase or real part.
- Current arrows show the **real vector phasor at 0°** of the selected harmonic,
  not a streamline calculation. Arrow length is scaled relative to the maximum.
- Phase has little meaning when the corresponding amplitude is near numerical
  noise. Inspect amplitudes and refinement before interpreting phase patterns.
- Click the field to select a temperature probe. Probe histories use the nearest
  node, rather than the cell-averaged display value.

## Export formats

Export model: JSON configuration. Export results: configuration, mesh coordinates,
full histories, terminal quantities and complex harmonic arrays. Spectrum CSV
contains terminal-voltage DC/1ω/2ω/3ω amplitude, phase, real and imaginary values.
Model edits disable result exports until a new run, to avoid associating old results
with changed controls. Results already displayed remain the previous computed run.

Indexing:

```
node index = j*(nx+1) + i        i=0..nx, j=0..ny
cell index = j*nx + i           i=0..nx-1, j=0..ny-1
```

j=0 is the bottom. Temperature/voltage histories use `[sample][node]`; Jx/Jy/qx/qy
use `[sample][cell]`. Harmonics use `[order][node or cell]` entries `{re,im}`.
Terminal quantities have one value per sample/order. The steady output contains
one spatial snapshot. `energyResidual` is heat leaving minus electrical power
absorbed, in W, and is reported for steady states only. `freeResidualWatts` is a
local semi-discrete thermal balance diagnostic, not a temperature-error estimate.

## Source files and editing

| File | Role |
|---|---|
| `src/index.template.html` | Five-tab layout |
| `src/styles.css` | Responsive styles |
| `src/app.js` | Controls, painting, canvas maps, plots, exports |
| `src/worker.js` | Background execution |
| `src/core/materials.js` | Material laws, adapted from the 1D version |
| `src/core/mesh2d.js` | Nodal grid, material cells, conservative links, boundary measures |
| `src/core/sparse.js` | Matrix-free graph solves with Dirichlet elimination |
| `src/core/solver2d.js` | 2D electric transport, thermal balance, nonlinear and periodic solves |
| `src/core/fourier.js` | Harmonic projection, using the same convention as the 1D version |
| `src/core/config2d.js` | Serializable configuration and model construction |
| `build.cjs` | Embeds source into the standalone index.html |

Edit the full files in `src/`, then rebuild with Node.js:

```bash
node build.cjs
```

No npm install is required for rebuilding or numerical tests.
The JS API is also usable from Node:

```js
const TE = require('./tests/load.cjs');
const config = TE.default2D();
config.mode = 'steady';
config.electrical.value = {bias: .1, amplitude: 0};
config.electrical.sourceRange = [.25, .75];
const result = TE.run2D(config);
console.log(result.terminalVoltage, result.energyResidual);
```

For custom numerical tolerances, construct `TE.from2DConfig(config)` and call
`solvePeriodic(f, {samples, maxPeriods, periodicAtol, periodicRtol, tolerance})`
or `solveSteady({tolerance, maxIterations, relaxation})` directly. Do not add
unrecognized numerical keys to GUI configuration expecting them to be used.

## Validation and remaining limits

Run:

```bash
node --test tests/solver.test.cjs tests/worker.test.cjs
```

24 numerical/packaging/worker tests pass. They include analytical 2D solutions,
current spreading, circulating currents in open circuit, interface Peltier,
Thomson, energy balances, mesh/time refinement, and agreement with the corrected
1D solver. See `VALIDATION.md` and `validation.txt`.

Full visual/click browser testing was not executed here because no browser binary
was available. The exact embedded worker source was tested with Node workers,
and an optional Playwright browser test is supplied. This does not substitute
for checking your browser's local-file/Blob policy or your target screen layout.

This model remains a rectangular, constant-depth, isotropic, perfect-contact
approximation. It does not support holes, curved or unstructured meshes, variable
depth, anisotropic tensors, contact resistances, front/back surface losses,
radiation or dynamic external electrical circuits. Material data and experimental
agreement must be validated for the intended application.
