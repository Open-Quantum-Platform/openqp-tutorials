# MRSF-TDDFT excited states and gradients

Most excited-state methods stumble exactly where photochemistry gets
interesting: near **conical intersections**, and for states with **open-shell or
diradical character**. Ordinary linear-response TDDFT can give qualitatively
wrong topology at an intersection and can badly imbalance singlet versus triplet
states. **MRSF-TDDFT** (Mixed-Reference Spin-Flip TDDFT) is built to fix this, and
it is OpenQP's workhorse excited-state method. You reach for it when you need
reliable excitation energies *and* analytic gradients — for excited-state
optimizations, minimum-energy conical-intersection searches, or as the electronic
engine behind nonadiabatic dynamics. This tutorial computes the three lowest MRSF
states of water (S0, S1, S2) and the analytic gradient of the first excited
state.

## A little theory

MRSF-TDDFT starts not from the closed-shell ground state but from a **high-spin
(triplet) ROHF reference**. From that reference it applies **spin-flip**
excitations — flipping one electron's spin while exciting it — to reach the target
singlet and triplet manifold. Plain spin-flip TDDFT does this from a *single*
reference determinant and suffers from spin contamination (the response states are
not clean spin eigenstates). MRSF's key move is to **mix two reference
determinants** (the M_s = +1 and M_s = -1 components of the triplet), which cancels
that contamination and yields balanced, well-defined singlet and triplet states —
precisely what you need when two states are about to cross.

Practically, three ingredients always travel together: a **triplet ROHF
reference**, the **spin-flip response** built on it, and a **half-and-half
functional** such as BHHLYP, which is the standard, well-benchmarked choice for
MRSF. In a `.oqp` deck the first two are carried by the single token `mrsf`, and
the third is the functional component of the route. The manifold that comes back
starts with the S0-like state, so you refer to the states simply as `S0`, `S1`,
`S2`. For the full derivation and benchmarks, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_mrsf.oqp`](inputs/h2o_mrsf.oqp) — water in the
6-31G\* basis, BHHLYP functional, three MRSF roots, plus the analytic gradient of
S1. Annotated:

```text
mrsf(nstate=3)/bhhlyp/6-31g* grad(S1)    # model(options)/functional/basis + gradient of S1
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Key points, line by line:

- **`mrsf(nstate=3)/bhhlyp/6-31g*`** is the whole electronic model. `mrsf` carries
  the high-spin **triplet ROHF reference** that MRSF spin-flips out of — you do not
  restate it, and you cannot accidentally pair MRSF with the wrong reference.
  `nstate=3` solves three roots; because the S0-like root is one of them, that
  gives S0, S1, and S2. `bhhlyp` is the half-and-half functional MRSF is normally
  run with; drop that component (`mrsf-tdhf/6-31g*`) for MRSF-TDHF.
- **`grad(S1)`** is the *driver*: energy **and** the analytic gradient, of the
  state named **physically**. You never work out that S1 is internal response root
  2 — that mapping is the format's job. Use `energy(S0)` for energies only, or
  `opt(S1)` to optimize the S1 minimum.
- **`geom`** holds the molecule inline; `geom="water.xyz"` reads it from a file
  instead.

Everything else — `method=tdhf`, `[tdhf] type=mrsf`, `[scf] type=rohf`,
`[scf] multiplicity=3`, `[properties] grad=2`, the Huckel initial guess, the SCF
thresholds — is implied by that one line or is a default, and is filled in for
you. Exact section calls such as `scf(conv=1e-8)` are the escape hatch when a
default is not what you want.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_mrsf.py`](inputs/h2o_mrsf.py). A single `job.theory.mrsf(...)` call
sets up the whole MRSF stack — under the hood it writes `[input] method=tdhf`,
`[scf] type=rohf multiplicity=3`, and `[tdhf] type=mrsf nstate=...` — and
`job.workflow.gradient(state=...)` selects `runtype=grad` and fills `[properties]
grad`:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf", silent=1)

# Water geometry (Angstrom), charge 0. The default multiplicity is fine here;
# job.theory.mrsf() sets the triplet ROHF *reference* multiplicity itself.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
)

# MRSF-TDDFT theory: BHHLYP/6-31G*, 3 response roots.
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)

# Gradient of state 2 (S1). state=... maps to [properties] grad and selects
# runtype=grad. For MRSF, state 1 is the lowest (S0-like) root, state 2 is S1.
job.workflow.gradient(state=2)

mol = job.run()

results = mol.get_results()
print("Reference (ROHF) energy:", results["energy"])
print("MRSF state total energies (S0, S1, S2):", results["td_energies"])
print("Gradient of S1 (Hartree/Bohr):", mol.get_grad())
```

Two conveniences worth noting: you do **not** set `multiplicity=3` on
`job.molecule(...)` — the triplet is the *reference* multiplicity and
`job.theory.mrsf()` sets it internally, so the molecule keeps its physical
(default) spin. And `job.workflow.gradient(state=2)` is what makes this a gradient
run rather than an energy-only run.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd mrsf-tddft/inputs
openqp h2o_mrsf.oqp
```

Python style:

```bash
cd mrsf-tddft/inputs
python h2o_mrsf.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

An MRSF gradient run reports three things: the **reference (ROHF) energy**, the
**total energy of each MRSF root**, and the **gradient of the requested state**.

- In the **log file** (`<project>.log`) look for the block of MRSF state energies
  (S0, S1, S2), the excitation energies relative to S0, and the printed Cartesian
  gradient of state 2.
- From **Python**:
  - `results["td_energies"]` (from `mol.get_results()`) is the list of **MRSF total
    energies**, ordered `[S0, S1, S2]`. The **excitation energy** of S1 is
    `td_energies[1] - td_energies[0]`, and of S2 is `td_energies[2] -
    td_energies[0]`.
  - `results["energy"]` is the reference-level energy printed alongside.
  - `mol.get_grad()` returns the **analytic gradient of S1** (state 2), in
    Hartree/Bohr, as one (x, y, z) triple per atom. This is the force you would feed
    to a geometry optimizer to relax the S1 surface.

To get *more* states, raise `nstate` in the route — `mrsf(nstate=5)/...` — or in
`job.theory.mrsf(...)`. To take the gradient of a different state, name it in the
driver: `grad(S2)`, or `job.workflow.gradient(state=...)` with that root's 1-based
index. For energies only, use `energy(S0)` / omit the
`job.workflow.gradient(...)` call.

## Manual

- MRSF-TDDFT workflow (reference, spin-flip response, functional choice, states):
  <https://open-quantum-platform.github.io/openqp-docs/>
- `[tdhf]` keyword reference (`type=mrsf`, `nstate`) and `[properties]` (`grad`):
  <https://open-quantum-platform.github.io/openqp-docs/>
