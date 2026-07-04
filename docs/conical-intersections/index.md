# Conical intersections: optimizing a MECI (MECP) with MRSF-TDDFT

When a molecule is electronically excited, the neat picture of nuclei gliding on a
single Born-Oppenheimer surface breaks down wherever two surfaces *touch*. Those
touching points are **conical intersections** — funnels through which an excited
molecule decays back to the ground state in femtoseconds, without emitting light.
They govern photostability, vision, photochemical switches, and the branching of
essentially every non-radiative photoreaction. This tutorial finds the lowest such
funnel between the ground state (S₀) and the first excited state (S₁) of
**ethylene**: its **minimum-energy conical intersection (MECI)**. For states of the
same spin the crossing is a genuine conical intersection; the identical machinery
finds a **minimum-energy crossing point (MECP)** between states of *different* spin
(e.g. a singlet/triplet seam), which is why the two names travel together.

## A little theory

Two adiabatic states become degenerate on a **seam** — a subspace of nuclear
configurations, not a single point. A MECI is the *lowest-energy* geometry on that
seam: the point you would actually reach after funnelling down S₁. Optimizing onto
it is harder than an ordinary minimization, because you must **minimize the energy
while simultaneously driving the S₀-S₁ gap to zero** — two competing objectives.

OpenQP's native optimizer offers the **penalty-function** method for this. Instead
of needing the (often ill-defined or expensive) non-adiabatic coupling and gradient-
difference vectors that span the branching plane, it minimizes a single smooth
objective built from the average energy of the two states plus a penalty that grows
as their gap widens. The penalty strength is ramped up over the optimization
(`pen_sigma`, `pen_incre`) until the gap closes to a tight tolerance (`energy_gap`)
at a stationary geometry. It is robust and needs only state energies and gradients.

The excited states themselves come from **MRSF-TDDFT** (Mixed-Reference Spin-Flip
TDDFT), which spin-flips out of a high-spin (triplet) ROHF reference and mixes two
reference determinants. That gives *balanced* descriptions of S₀ and S₁ right where
they cross — including the multiconfigurational, diradical character that plain
TDDFT gets wrong at a conical intersection. Ethylene is the textbook case: twisting
about the C=C bond and pyramidalizing one CH₂ carries the molecule to a twisted-
pyramidal S₀/S₁ funnel. For the underlying methods and the full input contract, see
the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/c2h4_mrsf_meci.inp`](inputs/c2h4_mrsf_meci.inp) —
twisted ethylene, BHHLYP/6-31G\*, MRSF-TDDFT, penalty-function MECI search.
Annotated:

```ini
[input]
system=
 C  -1.6699351346837055   0.1537249235528157  -1.5459803491111643   # twisted C=C
 C  -1.8079415266835852  -0.0386075716896284  -0.1602069788110266
 H  -2.6609567768367581   0.2572290722092156  -2.0290359598415040
 H  -1.2898503996116444  -0.7568524635289917  -2.0470428696820342
 H  -1.3096398768036397   0.6557118321425524   0.5396052278505126
 H  -2.3820842951209360  -0.7983813277099963   0.4308517619153288
runtype=meci            # optimize onto a conical-intersection seam
functional=bhhlyp       # half-and-half functional, standard for MRSF
charge=0
method=tdhf             # excited states via the TDHF/TDDFT engine (MRSF lives here)
basis=6-31g*

[scf]
type=rohf               # MRSF is built on a restricted-open-shell reference
maxit=30
multiplicity=3          # ...specifically a TRIPLET reference to spin-flip from

[tdhf]
type=mrsf               # Mixed-Reference Spin-Flip TDDFT
nstate=5                # solve 5 MRSF roots (S0 = root 1, S1 = root 2, ...)
maxit=30

[optimize]
lib=oqp                 # use OpenQP's native geometry optimizer
istate=1                # first  state of the crossing pair: S0
jstate=2                # second state of the crossing pair: S1
meci_search=penalty     # penalty-function branching-plane search
pen_sigma=2.0           # initial penalty strength
pen_incre=1.2           # factor by which the penalty is ramped each cycle
energy_gap=2e-3         # target S0-S1 gap (Hartree) -> defines "degenerate"
rmsd_grad=2e-3          # convergence: RMS gradient
max_grad=4e-3           # convergence: max gradient component
rmsd_step=4e-3          # convergence: RMS step
max_step=8e-3           # convergence: max step component
maxit=30                # max optimization cycles

[oqp]
coordsys=auto           # let the optimizer choose the coordinate system
trust=0.15              # trust-radius (Bohr/rad) for each optimization step
```

Key points:

- **`runtype=meci`** is what turns this into a crossing-point optimization rather
  than an ordinary geometry optimization or single-point. The same optimizer also
  drives `optimize`, `ts`, `neb`, `irc`, and related run types; `meci` selects the
  crossing-seam search.
- **The excited-state block is the same MRSF recipe used everywhere else:**
  `[scf] type=rohf multiplicity=3` sets up the triplet reference, and
  `[tdhf] type=mrsf nstate=5` spin-flips to the balanced singlet/triplet manifold.
  You need `nstate` large enough to contain both states of interest — here
  `nstate=5` comfortably covers roots 1 and 2.
- **`istate` / `jstate` name the two states to bring together** (1-based over the
  MRSF roots). `istate=1, jstate=2` is the S₀/S₁ pair. To find an S₀/S₂ funnel you
  would set `jstate=3`; for a spin-different MECP you would pick roots of different
  spin from the manifold.
- **`meci_search=penalty`** picks the penalty-function algorithm described above.
  `pen_sigma` is its starting strength and `pen_incre` the per-cycle ramp factor;
  `energy_gap` is the gap (in Hartree) below which the two states count as
  degenerate. Loosening `energy_gap` gives a quicker but looser crossing point.
- **The four `*_grad` / `*_step` thresholds plus `energy_gap` together define
  convergence:** the optimizer stops when the geometry is stationary *and* the gap
  is closed. `maxit` caps the cycles.
- **`[oqp]`** configures the native optimizer engine itself: `coordsys=auto` lets
  it choose internal vs. Cartesian coordinates, and `trust=0.15` sets the step-size
  trust radius.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/c2h4_mrsf_meci.py`](inputs/c2h4_mrsf_meci.py). `job.theory.mrsf(...)` fills
the `[scf]`/`[tdhf]` blocks (the triplet ROHF reference and the MRSF roots), and
`job.workflow.meci(...)` sets `runtype=meci` and fills the `[optimize]` and `[oqp]`
sections in one call — the `coordsys`/`trust` arguments are routed to `[oqp]`, the
rest to `[optimize]`.

```python
from oqp.openqp import OpenQP

geometry = """
C  -1.6699351346837055   0.1537249235528157  -1.5459803491111643
C  -1.8079415266835852  -0.0386075716896284  -0.1602069788110266
H  -2.6609567768367581   0.2572290722092156  -2.0290359598415040
H  -1.2898503996116444  -0.7568524635289917  -2.0470428696820342
H  -1.3096398768036397   0.6557118321425524   0.5396052278505126
H  -2.3820842951209360  -0.7983813277099963   0.4308517619153288
"""

job = OpenQP("c2h4_mrsf_meci", silent=1)

# Neutral molecule; MRSF builds its states on a triplet ROHF reference.
job.molecule(geometry, charge=0)

# MRSF-TDDFT on a BHHLYP/6-31G* reference, solving 5 roots.
# Root 1 is the ground state S0, root 2 is S1.
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=5)

# Crossing-point search between state 1 (S0) and state 2 (S1).
# lib="oqp" selects the native optimizer; coordsys/trust go to [oqp];
# energy_gap is the S0-S1 degeneracy target that makes this an MECI.
job.workflow.meci(
    lib="oqp",
    istate=1,
    jstate=2,
    meci_search="penalty",
    pen_sigma=2.0,
    pen_incre=1.2,
    energy_gap=2e-3,
    rmsd_grad=2e-3,
    max_grad=4e-3,
    rmsd_step=4e-3,
    max_step=8e-3,
    maxit=30,
    coordsys="auto",
    trust=0.15,
)

mol = job.run()

results = mol.get_results()
print("MRSF state energies at the optimized geometry:", results["td_energies"])
```

Every argument maps one-to-one to a keyword in the `.inp` deck, so the two scripts
run the *same* optimization and land on the same MECI geometry.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd conical-intersections/inputs
openqp c2h4_mrsf_meci.inp
```

Python style:

```bash
cd conical-intersections/inputs
python c2h4_mrsf_meci.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same result.

## Reading the output

The optimization prints one line per cycle; the numbers to watch are the
**energies of states `istate` and `jstate` and the gap between them**, alongside the
gradient/step norms. The run has succeeded when, on the final cycle:

- the **S₀-S₁ energy gap** has dropped below `energy_gap` (2e-3 Ha ≈ 0.05 eV) — the
  two states are degenerate, and
- the **gradient and step** norms are below their `*_grad` / `*_step` thresholds —
  the geometry is stationary on the seam.

That converged geometry *is* the MECI: the twisted-pyramidal ethylene funnel.

- In the **log file** (`<project>.log`) look for the per-cycle state energies, the
  gap, and the convergence flags; the final block reports the optimized MECI
  geometry and the two (now nearly equal) state energies.
- From **Python**, `mol.get_results()["td_energies"]` returns the MRSF state energies
  at the optimized geometry — the tutorial script prints exactly this. The first two
  entries (S₀ and S₁) should agree to within `energy_gap`, confirming the crossing.
  `mol.get_scf_energy()` gives the underlying ROHF reference energy if you need it.

To explore further: raise `jstate` to hunt a higher crossing (S₀/S₂), loosen or
tighten `energy_gap` to trade speed for how closed the seam must be, or start from a
different guess geometry to locate a different funnel on the same seam.

## Manual

- Geometry-optimization / crossing-point workflow (MECI/MECP, `istate`/`jstate`,
  the native `lib=oqp` optimizer):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/optimization/>
- `[optimize]` keyword reference (`meci_search`, `pen_sigma`, `pen_incre`,
  `energy_gap`, the convergence thresholds):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/optimize/>
- `[oqp]` keyword reference (`coordsys`, `trust` — configuring the native optimizer):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/oqp/>
- MRSF-TDDFT excited states (`[scf]` triplet reference, `[tdhf] type=mrsf`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/tdhf/>
```
