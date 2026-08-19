# Geometry optimization and transition states: finding stationary points

A single-point energy tells you what a molecule costs at *one* fixed geometry.
Most chemistry questions, though, are about **stationary points** on the
potential energy surface (PES) — the shapes a molecule actually adopts. A
**geometry optimization** walks downhill to the nearest **minimum** (an
equilibrium structure: a reactant, product, or stable conformer), while a
**transition-state (TS) search** climbs to the first-order **saddle point** that
separates two minima — the bottleneck geometry whose energy sets a reaction
barrier. You run these before almost anything else: frequencies, reaction
energies, and barrier heights are all defined *at* optimized geometries. This
tutorial does both on small molecules — a minimum for water, a saddle point for
the HCN → HNC isomerization.

## A little theory

Both jobs are driven by the energy gradient **g = ∂E/∂x** (the force) and,
implicitly or explicitly, the Hessian **H = ∂²E/∂x²** (its curvature). A
**minimization** takes quasi-Newton steps `Δx = -H⁻¹g` that lower the energy
until the gradient vanishes and the Hessian is positive-definite — every
direction curves upward. A **TS search** looks for the same `g = 0` condition
but at a point with **exactly one negative Hessian eigenvalue**: downhill along
the reaction coordinate, uphill in every other direction. That one extra
constraint is what makes a saddle point harder to find than a minimum, so TS
searches lean more on a good starting guess and on eigenvector-following steps
that deliberately go *up* along the reaction mode.

Two practical knobs recur below. The **coordinate system** the optimizer steps
in — Cartesian, delocalized internal coordinates (**DLC**), or translation-
rotation internal coordinates (**TRIC**) — strongly affects how fast it
converges; internal-coordinate systems remove the redundant translations and
rotations and follow chemical bonds, so they usually beat raw Cartesians. The
**trust radius** caps how far a single step may move, keeping the quadratic model
honest. OpenQP ships two optimizer backends: the native `oqp` optimizer, which the
concise `.oqp` format always uses, and the external
[geomeTRIC](https://geometric.readthedocs.io/) library, still selectable from the
Python API and from legacy `.inp` decks. For the full contract see the
[optimization workflow page](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

### Minimum: water (`inputs/h2o_optimize.oqp`)

The first deck relaxes water to its nearest minimum. Annotated:

```text
rks/bhhlyp/6-31g* opt            # BHHLYP Kohn-Sham gradients; relax to the nearest MINIMUM
geom="""
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
"""
```

Line by line:

- **`rks/bhhlyp/6-31g*`** is the level of theory the gradient is taken at: a
  restricted Kohn-Sham reference with the BHHLYP functional. Drop the middle
  component (`rhf/6-31g*`) for a Hartree-Fock optimization.
- **`opt(...)` is the driver**, and it is what turns a single point into a
  **minimization**. With a ground-state model there is no state to name; on a
  response route you would write `opt(S1)` to relax on an excited surface.
- **Nothing else is needed.** The concise format always uses the **native OpenQP
  optimizer**, so there is no backend to select and no separate backend section.
  Its step controls are options of the driver and default to sensible values —
  an automatic internal-coordinate choice (`coordsys=auto`, which starts in DLC),
  a trust radius (maximum step) of 0.2, and at most 30 cycles. Write one only to
  change it, for example `opt(coordsys=tric,trust=0.1,maxit=60)`.

### Transition state: HCN → HNC (`inputs/hcn_ts.oqp`)

The second deck searches for the **saddle point** of the HCN ⇌ HNC
isomerization. The only structural difference from a minimization is the driver
name.

```text
rks/bhhlyp/3-21g ts(S0,maxit=50,coordsys=dlc,trust=0.05)  # small basis; first-order SADDLE-POINT search
geom="""
C   0.0000000000   0.0000000000   0.0000000000
N   0.0000000000   0.0000000000   1.1700000000
H  -1.1000000000   0.0000000000   0.0000000000
"""
```

What changes relative to the water minimization:

- **`ts(...)` replaces `opt(...)`**, switching the driver from downhill
  minimization to a first-order saddle-point search (native P-RFO). That single
  token is the whole difference in intent.
- **`coordsys=dlc`** works in delocalized internal coordinates, and a deliberately
  small **`trust=0.05`** keeps each step inside the quadratic region — TS steps are
  easy to overshoot.
- **`maxit=50`** — saddle points usually take more cycles than minima, so the
  ceiling is raised above the default 30.

`ts(...)` also takes `follow` (which mode to follow uphill) and `hessian`
(how to seed the initial Hessian); see the manual for the full list.

> **Note.** Earlier versions of this deck ran the saddle-point search through the
> external **geomeTRIC** backend, with its own `tmax`, `convergence_set=GAU`, and
> `hessian=never` vocabulary. The concise `.oqp` format deliberately hides
> optimizer-backend selection and always uses the native engine, so those
> geomeTRIC-only keys have no `.oqp` spelling. `opt(...)` — but not `ts(...)` —
> still accepts an explicit `lib=geometric`. If you need geomeTRIC for a TS
> search, use a legacy `.inp` deck with `[optimize] lib=geometric`.

## Python style

Each deck has a one-to-one Python twin using the compact `OpenQP` scripting
interface. `job.theory.dft(...)` sets the reference, and `job.workflow.optimize(...)`
/ `job.workflow.ts(...)` fill the `[optimize]` section and route the extra
keywords to the right backend section automatically.

### Minimum: `inputs/h2o_optimize.py`

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_optimize", silent=1)

# Same starting geometry as the .oqp (Angstrom); charge 0, closed-shell singlet.
job.molecule(
    """
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
""",
    charge=0,
    multiplicity=1,
)

# BHHLYP Kohn-Sham reference with the 6-31g* basis, matching the deck's route.
job.theory.dft(functional="bhhlyp", basis="6-31g*", reference="rhf")

# Geometry optimization on the native optimizer.
#   lib="oqp"      -> [optimize] lib=oqp   (native backend)
#   istate=0       -> [optimize] istate=0  (ground state)
#   maxit=30       -> [optimize] maxit=30
#   coordsys/trust -> routed to the [oqp] backend section automatically
job.workflow.optimize(
    lib="oqp",
    istate=0,
    maxit=30,
    coordsys="tric",
    trust=0.2,
)

mol = job.run()

# The final SCF energy is the energy at the optimized minimum.
print("Optimized SCF energy:", mol.get_scf_energy())
print("Optimized geometry (Bohr):", mol.get_system())
print(mol.get_results())
```

The key mapping: `job.workflow.optimize(...)` is `runtype=optimize`; `lib`,
`istate`, and `maxit` fill `[optimize]`; and because `lib="oqp"`, `coordsys` and
`trust` are dispatched to the `[oqp]` section for you.

### Transition state: `inputs/hcn_ts.py`

The TS script is identical in shape — the single change that flips minimization
into a saddle-point search is calling **`job.workflow.ts(...)`** (which selects
`runtype=ts`) with `lib="geometric"`, so the extra keywords land in `[geometric]`:

```python
from oqp.openqp import OpenQP

job = OpenQP("hcn_ts", silent=1)

# A bent starting guess near the HCN <-> HNC saddle point (Angstrom).
job.molecule(
    """
C   0.0000000000   0.0000000000   0.0000000000
N   0.0000000000   0.0000000000   1.1700000000
H  -1.1000000000   0.0000000000   0.0000000000
""",
    charge=0,
    multiplicity=1,
)

# BHHLYP Kohn-Sham reference with a small 3-21g basis (fast).
job.theory.dft(functional="bhhlyp", basis="3-21g", reference="rhf")

# Transition-state search (runtype=ts) on the geomeTRIC backend.
#   lib="geometric"          -> [optimize] lib=geometric
#   istate=0                 -> [optimize] istate=0 (ground state)
#   coordsys/trust/tmax/...  -> routed to the [geometric] backend section
job.workflow.ts(
    lib="geometric",
    istate=0,
    maxit=50,
    coordsys="dlc",
    trust=0.05,
    tmax=0.1,
    hessian="never",
    convergence_set="GAU",
    prefix="hcn_ts",
)

mol = job.run()

# Energy at the converged saddle point.
print("TS SCF energy:", mol.get_scf_energy())
print(mol.get_results())
```

## Run it

Both styles produce the same result; run from the `inputs/` folder.

Water minimization:

```bash
cd geometry-optimization/inputs
openqp h2o_optimize.oqp       # input-file style
python h2o_optimize.py        # Python-API style
```

HCN → HNC transition state:

```bash
cd geometry-optimization/inputs
openqp hcn_ts.oqp             # input-file style
python hcn_ts.py              # Python-API style
```

Both need OpenQP installed (`pip install openqp`). The `.oqp` decks run entirely
on the native optimizer; the TS *Python* script above selects the external
**geomeTRIC** backend, which needs `pip install geometric`.

## Reading the output

An optimization run prints one line per geometry cycle — the energy, the maximum
and RMS gradient, and the step taken — and then reports convergence. The numbers
that matter:

- **`mol.get_scf_energy()`** is the SCF (here BHHLYP Kohn-Sham) energy **at the final
  geometry** — the energy of the optimized minimum, or the energy at the
  converged saddle point. This is the number you carry forward into barrier
  heights and reaction energies.
- **`mol.get_system()`** returns the **optimized geometry** (Cartesian
  coordinates, in Bohr) — the relaxed structure itself.
- **`mol.get_results()`** is the full results dictionary (energy, gradient, and
  run metadata), also written to `<project>.json`.
- In the **log file** (`<project>.log`) watch the gradient columns fall toward
  the convergence thresholds and look for the "converged" banner. If it hits
  `maxit` first, the geometry did **not** converge — tighten the starting guess
  or adjust the trust radius.

Two sanity checks specific to each job type. For the **minimum**, the water
geometry should relax to its equilibrium bond length and H–O–H angle with the
gradient driven to zero. For the **transition state**, a genuine first-order
saddle point has **exactly one imaginary vibrational frequency** — confirm it by
running a frequency (Hessian) calculation at the converged TS geometry; the
single imaginary mode should visibly connect HCN to HNC (the H migrating across
the C≡N axis).

## Manual

- Optimization / transition-state workflow (backends, coordinate systems, the
  `[optimize]` contract): <https://open-quantum-platform.github.io/openqp-docs/>
- `[optimize]` keyword reference (`lib`, `istate`, `maxit`) and the backend
  sections `[oqp]` / `[geometric]` (`coordsys`, `trust`, `tmax`, `hessian`,
  `convergence_set`): <https://open-quantum-platform.github.io/openqp-docs/>
- Running OpenQP from Python (the `job.workflow.optimize(...)` /
  `job.workflow.ts(...)` idiom):
  <https://open-quantum-platform.github.io/openqp-docs/>
