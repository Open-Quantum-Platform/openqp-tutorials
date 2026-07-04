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
the HCN → HNC isomerization — and shows OpenQP's two optimizer backends.

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
honest. OpenQP ships two optimizer backends you can select per job: the native
`oqp` optimizer and the external [geomeTRIC](https://geometric.readthedocs.io/)
library. For the full contract see the
[optimization workflow page](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

### Minimum: water (`inputs/h2o_optimize.inp`)

The first deck relaxes water to its nearest minimum with the **native `oqp`**
optimizer. Annotated:

```ini
[input]
system=
   O  -0.0000000000   0.0000000000  -0.0410615540   # starting geometry (Angstrom)
   H  -0.5331943294   0.5331943294  -0.6144692230
   H   0.5331943294  -0.5331943294  -0.6144692230
charge=0
functional=bhhlyp       # functional is set but ignored — method=hf uses HF exchange only
basis=6-31g*
runtype=optimize        # relax to the nearest MINIMUM
method=hf               # Hartree-Fock energy and gradient

[scf]
type=rhf                # closed-shell restricted HF reference
multiplicity=1          # singlet

[optimize]
lib=oqp                 # optimizer backend: native OpenQP
istate=0                # state to optimize: 0 = HF/DFT ground state
maxit=30                # max optimization cycles

[oqp]                   # backend section for lib=oqp
coordsys=tric           # step in translation-rotation internal coordinates
trust=0.2               # trust radius (max step size)
```

Section by section:

- **`[input]`** holds the starting geometry and the level of theory. `runtype=optimize`
  is what turns a single point into a **minimization**; `method=hf` picks a
  Hartree-Fock energy/gradient. (`functional=bhhlyp` is present but inert here —
  with `method=hf` OpenQP uses exact HF exchange and no DFT correlation.)
- **`[scf]`** defines the reference the gradient is taken on: `type=rhf`,
  closed-shell `multiplicity=1`.
- **`[optimize]`** is the backend-independent driver. **`lib=oqp`** selects the
  native optimizer; **`istate=0`** says optimize the ground state (raise it to
  follow an excited state); **`maxit=30`** caps the number of geometry cycles.
- **`[oqp]`** carries the settings *specific to* the `lib=oqp` backend:
  **`coordsys=tric`** (TRIC internal coordinates) and **`trust=0.2`** (the step
  cap, in the optimizer's internal units).

### Transition state: HCN → HNC (`inputs/hcn_ts.inp`)

The second deck searches for the **saddle point** of the HCN ⇌ HNC
isomerization, using the **geomeTRIC** backend. The only structural difference
from a minimization is `runtype=ts` and a different backend section.

```ini
[input]
system=
   C   0.0000000000   0.0000000000   0.0000000000   # bent guess near the saddle (Angstrom)
   N   0.0000000000   0.0000000000   1.1700000000
   H  -1.1000000000   0.0000000000   0.0000000000
charge=0
functional=bhhlyp       # inert with method=hf (see above)
basis=3-21g             # small, fast basis for the demo
runtype=ts              # search for a first-order SADDLE POINT
method=hf

[scf]
type=rhf
multiplicity=1
maxit=30                # max SCF iterations per point

[optimize]
lib=geometric           # optimizer backend: external geomeTRIC
istate=0                # ground-state saddle point
maxit=50                # TS searches often need more cycles than minima

[geometric]             # backend section for lib=geometric
coordsys=dlc            # delocalized internal coordinates
trust=0.05              # small trust radius — TS steps must stay in the quadratic region
tmax=0.1                # hard cap on the maximum single step
hessian=never           # do not build/update an exact Hessian (use the model Hessian)
convergence_set=GAU     # Gaussian-style convergence thresholds
prefix=hcn_ts           # base name for geomeTRIC's own output files
```

What changes relative to the water minimization:

- **`runtype=ts`** switches the driver from downhill minimization to a
  first-order saddle-point search.
- **`[optimize] lib=geometric`** selects the external geomeTRIC engine, so the
  backend-specific keys live in a **`[geometric]`** section instead of `[oqp]`.
- The geomeTRIC keys are its own vocabulary: **`coordsys=dlc`** (delocalized
  internals), a deliberately small **`trust=0.05`** with a hard **`tmax=0.1`**
  step cap (TS steps are easy to overshoot), **`hessian=never`** (rely on the
  cheap model Hessian rather than computing an exact one), **`convergence_set=GAU`**
  (the familiar Gaussian force/displacement thresholds), and **`prefix=hcn_ts`**
  for the names of geomeTRIC's log files.
- **`maxit=50`** in `[optimize]` — saddle points usually take more cycles than
  minima, so the ceiling is raised.

The pattern is the same for both jobs: `[optimize]` chooses the backend and the
state, and a backend-named section (`[oqp]` **or** `[geometric]`) carries that
backend's step-control keywords.

## Python style

Each deck has a one-to-one Python twin using the compact `OpenQP` scripting
interface. `job.theory.hf(...)` sets the reference, and `job.workflow.optimize(...)`
/ `job.workflow.ts(...)` fill the `[optimize]` section and route the extra
keywords to the right backend section automatically.

### Minimum: `inputs/h2o_optimize.py`

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_optimize", silent=1)

# Same starting geometry as the .inp (Angstrom); charge 0, closed-shell singlet.
job.molecule(
    """
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
""",
    charge=0,
    multiplicity=1,
)

# Hartree-Fock reference (method=hf) with the 6-31g* basis.
job.theory.hf(basis="6-31g*")

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

# Hartree-Fock reference with a small 3-21g basis (fast).
job.theory.hf(basis="3-21g")

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
openqp h2o_optimize.inp       # input-file style
python h2o_optimize.py        # Python-API style
```

HCN → HNC transition state:

```bash
cd geometry-optimization/inputs
openqp hcn_ts.inp             # input-file style
python hcn_ts.py              # Python-API style
```

Both need OpenQP installed (`pip install openqp`); the TS deck additionally uses
the external **geomeTRIC** backend (`pip install geometric`).

## Reading the output

An optimization run prints one line per geometry cycle — the energy, the maximum
and RMS gradient, and the step taken — and then reports convergence. The numbers
that matter:

- **`mol.get_scf_energy()`** is the SCF (here HF) energy **at the final
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
