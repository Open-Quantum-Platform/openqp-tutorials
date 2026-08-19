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

OpenQP's native optimizer handles this without the (often ill-defined or
expensive) non-adiabatic coupling and gradient-difference vectors that span the
branching plane: it minimizes a single smooth objective built from the average
energy of the two states plus a term that grows as their gap widens, until the
gap closes to a tight tolerance at a stationary geometry. For two states the
default is an **augmented Lagrangian** (`auglag`), which converges the gap
itself; the older **penalty-function** method (`algorithm=penalty`, with a
ramped penalty strength `sigma`/`pen_incre`) and the multistate **BaekA**
algorithm are available as well. All of them need only state energies and
gradients.

The excited states themselves come from **MRSF-TDDFT** (Mixed-Reference Spin-Flip
TDDFT), which spin-flips out of a high-spin (triplet) ROHF reference and mixes two
reference determinants. That gives *balanced* descriptions of S₀ and S₁ right where
they cross — including the multiconfigurational, diradical character that plain
TDDFT gets wrong at a conical intersection. Ethylene is the textbook case: twisting
about the C=C bond and pyramidalizing one CH₂ carries the molecule to a twisted-
pyramidal S₀/S₁ funnel. For the underlying methods and the full input contract, see
the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/c2h4_mrsf_meci.oqp`](inputs/c2h4_mrsf_meci.oqp) —
twisted ethylene, BHHLYP/6-31G\*, MRSF-TDDFT, native MECI search. Annotated:

```text
mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1)    # five MRSF roots; S0/S1 crossing search
geom="""
C   -1.6699351346837055    0.1537249235528157   -1.5459803491111643
C   -1.8079415266835852   -0.0386075716896284   -0.1602069788110266
H   -2.6609567768367581    0.2572290722092156   -2.0290359598415040
H   -1.2898503996116444   -0.7568524635289917   -2.0470428696820342
H   -1.3096398768036397    0.6557118321425524    0.5396052278505126
H   -2.3820842951209360   -0.7983813277099963    0.4308517619153288
"""
```

That one line is the whole request. Every algorithm and convergence control
has a default, so none is written here; the list below says what they are and
how to change one.

Key points:

- **`meci(...)` is the driver**, and it is what turns this into a crossing-point
  optimization rather than an ordinary geometry optimization or single point. The
  same native optimizer also backs `opt`, `ts`, `neb`, `irc`, and `mecp`; naming
  `meci` selects the crossing-seam search.
- **`S0,S1` name the two states to bring together — physically.** You write the
  states you mean; the format maps them to the internal MRSF roots. To find an
  S₀/S₂ funnel, write `meci(S0,S2,algorithm=penalty)` — the multistate BaekA
  algorithm requires *consecutive* roots in one spin manifold and says so
  rather than quietly optimizing the wrong pair. For a spin-different crossing,
  name states of different spin (`mecp(S0,T1)`) and use the `mecp` driver.
- **`mrsf(nstate=5)/bhhlyp/6-31g*`** is the same MRSF recipe used everywhere else:
  the triplet ROHF reference and the spin-flip to a balanced singlet/triplet
  manifold, all carried by the one `mrsf` token. `nstate` must be large enough to
  contain both states of interest — `nstate=5` comfortably covers S0 and S1.
- **The algorithm is chosen for you.** Two states default to the augmented
  Lagrangian; `meci(S0,S1,algorithm=penalty)` selects the penalty-function
  method instead, whose starting strength and per-cycle ramp are `sigma` and
  `pen_incre`. `gap` (default `1e-5` Hartree) is the energy gap below which the
  two states count as degenerate; loosening it, e.g. `gap=2e-3`, gives a quicker
  but looser crossing point. The crossing drivers use these short public names —
  `algorithm`, `sigma`, `alpha`, `gap` — rather than the legacy
  `meci_search`/`pen_sigma`/`energy_gap` spellings.
- **Convergence is defined by the `rmsd_grad`/`max_grad`/`rmsd_step`/`max_step`
  thresholds plus `gap`:** the optimizer stops when the geometry is stationary
  *and* the gap is closed. `maxit` (default 30) caps the cycles. Add any of them
  to the driver call only to loosen or tighten the search.
- **`coordsys` and `trust`** configure the native optimizer engine itself —
  internal or Cartesian coordinates (`auto` starts in DLC) and the step-size
  trust radius (0.2). In the concise format these are options of the driver,
  not a separate engine section.
- **SCF and response-solver controls** such as `scf(conv=1e-8)` or
  `tdhf(maxit=100)` are exact section calls; they are not needed here.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/c2h4_mrsf_meci.py`](inputs/c2h4_mrsf_meci.py). `job.theory.mrsf(...)` fills
the `[scf]`/`[tdhf]` blocks (the triplet ROHF reference and the MRSF roots), and
`job.workflow.meci(...)` sets `runtype=meci` and fills `[optimize]` (and, for
engine controls such as `coordsys`/`trust`, `[oqp]`) in one call.

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

# Crossing-point search between state 1 (S0) and state 2 (S1) with the native
# optimizer and its default (augmented-Lagrangian) algorithm and thresholds.
job.workflow.meci(istate=1, jstate=2)

mol = job.run()

results = mol.get_results()
print("MRSF state energies at the optimized geometry:", results["td_energies"])
```

The arguments map one-to-one to the `.oqp` deck, so the two scripts run the
*same* optimization and land on the same MECI geometry.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd conical-intersections/inputs
openqp c2h4_mrsf_meci.oqp
```

Python style:

```bash
cd conical-intersections/inputs
python c2h4_mrsf_meci.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same result.

## Reading the output

The optimization prints one line per cycle; the numbers to watch are the
**energies of the two states named in the driver and the gap between them**,
alongside the gradient/step norms. The run has succeeded when, on the final cycle:

- the **S₀-S₁ energy gap** has dropped below `gap` (default 1e-5 Ha) — the
  two states are degenerate, and
- the **gradient and step** norms are below their `*_grad` / `*_step` thresholds —
  the geometry is stationary on the seam.

That converged geometry *is* the MECI: the twisted-pyramidal ethylene funnel.

- In the **log file** (`<project>.log`) look for the per-cycle state energies, the
  gap, and the convergence flags; the final block reports the optimized MECI
  geometry and the two (now nearly equal) state energies.
- From **Python**, `mol.get_results()["td_energies"]` returns the MRSF state energies
  at the optimized geometry — the tutorial script prints exactly this. The first two
  entries (S₀ and S₁) should agree to within `gap`, confirming the crossing.
  `mol.get_scf_energy()` gives the underlying ROHF reference energy if you need it.

To explore further: name a higher pair to hunt a higher crossing
(`meci(S0,S2,algorithm=penalty)`), loosen or tighten `gap` to trade speed for how
closed the seam must be, or start from a different guess geometry to locate a
different funnel on the same seam.

## Manual

- Geometry-optimization / crossing-point workflow (MECI/MECP, `istate`/`jstate`,
  the native `lib=oqp` optimizer):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/optimization/>
- `[optimize]` keyword reference — the legacy names behind the driver's
  `algorithm` / `sigma` / `pen_incre` / `gap` and the convergence thresholds:
  <https://open-quantum-platform.github.io/openqp-docs/keywords/optimize/>
- `[oqp]` keyword reference (`coordsys`, `trust` — configuring the native optimizer):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/oqp/>
- MRSF-TDDFT excited states (`[scf]` triplet reference, `[tdhf] type=mrsf`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/tdhf/>
```
