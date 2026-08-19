# Spin-orbit coupling between MRSF states

Non-relativistic quantum chemistry treats spin as a good quantum number:
singlets and triplets live in separate worlds and never talk to each other. The
**spin-orbit coupling (SOC)** operator — a relativistic effect where an
electron's spin feels the magnetic field of its own orbital motion — is what
breaks that wall. Its matrix elements between singlet and triplet states set the
rate of **intersystem crossing** and **phosphorescence**, so you compute them
whenever spin-forbidden processes matter: photosensitizers, OLED emitters,
triplet harvesting, heavy-atom photochemistry. This tutorial computes the SOC
matrix elements between the **MRSF-TDDFT** singlet and triplet states of a water
molecule — a single-point calculation that is the same electronic-structure step
sitting inside a SOC-NAMD dynamics run, just without the nuclear motion.

## A little theory

**MRSF-TDDFT** (Mixed-Reference Spin-Flip TDDFT) reaches the excited states by
spin-flipping out of a **high-spin (triplet) ROHF reference**. Mixing two
reference determinants removes the spin contamination that plagues plain
spin-flip TDDFT and delivers balanced **singlet *and* triplet** roots from one
calculation — exactly the manifold you need when you want couplings *between*
the two spin blocks.

The coupling itself comes from the **Breit-Pauli spin-orbit operator**. OpenQP
can evaluate it at two levels: the one-electron term alone (`soc_2e=0`), or the
one-electron term plus a **mean-field two-electron** contribution (`soc_2e=1`),
which is much more accurate and still affordable. Because these integrals are
formulated over **Cartesian** Gaussians, the two-electron SOC path requires
`ispher=false`. Scalar-relativistic effects on the core are folded in separately
through the **Douglas-Kroll-Hess (DKH)** transformation (`scal_rel`). OpenQP
diagonalizes the spin-adiabatic Hamiltonian built from the spin-free energies
plus the SOC matrix, and reports the resulting eigenvalues in cm⁻¹. For the full
derivation and the operator definitions, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_soc.oqp`](inputs/h2o_soc.oqp) — water in the
6-31G(2df,p) basis, a triplet-ROHF MRSF reference, SOC over the singlet/triplet
manifold. Annotated:

```text
mrsf(nstate=12)/bhhlyp/6-31g(2df,p) ispher=false soc   # 12 MRSF roots, Cartesian AOs, SOC driver
scf(scal_rel=2,maxit=200)                               # DKH2 scalar relativity; raised iteration cap
geom="""
O    0.000000000000   0.000000000000   0.000000000000
H    0.772597940000   0.555677850000   0.000000000000
H   -0.773127700000   0.555677850000   0.000000000000
"""
```

Key points:

- **`soc(...)` is the driver.** It selects the spin-orbit-coupling workflow:
  converge the reference, solve the MRSF roots, build the SOC matrix, and
  diagonalize.
- **`mrsf(nstate=12)/bhhlyp/6-31g(2df,p)`** is the response model. `mrsf` carries
  the high-spin triplet ROHF reference MRSF spin-flips from — not a state you care
  about, but the mathematical starting point that produces balanced singlets and
  triplets. `nstate=12` requests twelve response roots; from these the SOC workflow
  forms the singlet and triplet states internally and couples them, so you do not
  enumerate the spin blocks yourself.
- **The mean-field two-electron SOC is on by default** (`soc_2e=1`) on top of
  the one-electron Breit-Pauli term. Because that path is built over Cartesian
  Gaussians, it **requires `ispher=false`**; write `soc(soc_2e=0)` if you only
  want the one-electron term. To size the singlet and triplet blocks independently, drop
  `nstate` from the route and give `soc(ns=...,nt=...)` instead — asking for both
  at once is rejected rather than silently resolved.
- **`scf(scal_rel=2)`** applies the second-order DKH scalar-relativistic correction
  to the core (`scal_rel=1` is first-order only; `0` turns it off), alongside a
  raised iteration cap. These are `[scf]` controls, so they travel in an exact
  section call rather than the route.
- **A basis name may contain parentheses and commas** — `6-31g(2df,p)` parses as
  one route component, because the route is split only at top level.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_soc.py`](inputs/h2o_soc.py). `job.theory.mrsf(...)` sets the triplet
ROHF reference and the MRSF roots; `job.workflow.soc(...)` selects
`runtype=soc` and fills in `soc_2e` and the `[scf] scal_rel` core treatment.

```python
from oqp.openqp import OpenQP

# Cartesian AOs (ispher=false): 2e SOC is not available with spherical AOs.
job = OpenQP("h2o_soc", silent=1, input={"ispher": "false"})

# Same water geometry as the input deck (Angstrom).
job.molecule(
    "8  0.000000  0.000000  0.000000; "
    "1  0.772598  0.555678  0.000000; "
    "1 -0.773128  0.555678  0.000000",
    charge=0,
)

# MRSF-TDDFT theory: triplet ROHF reference + [tdhf] type=mrsf, 12 roots.
job.theory.mrsf(functional="bhhlyp", basis="6-31G(2df,p)", nstate=12)

# SOC workflow: sets runtype=soc, soc_2e, and [scf] scal_rel.
#   soc_2e=1  -> 1e Breit-Pauli + mean-field 2e SOC
#   scal_rel=2 -> DKH scalar-relativistic core (helper default is 2)
job.workflow.soc(soc_2e=1, scal_rel=2)

mol = job.run()

# SOC eigenvalues (spin-adiabatic states) in cm^-1.
print("SOC eigenvalues (cm^-1):", mol.get_soc())

# Full JSON-friendly summary: td_singlet_energies, td_triplet_energies, soc, ...
print(mol.get_results())
```

A few things to notice:

- **`ispher="false"` is passed to the constructor** through `input={...}`, since
  it is an `[input]`-level flag and (as above) the two-electron SOC path needs
  Cartesian AOs.
- **`job.theory.mrsf(...)`** carries the same MRSF settings as the deck — the
  triplet ROHF reference is implied by the MRSF helper, and `nstate=12` requests
  the same twelve roots. The functional and basis match `bhhlyp` /
  `6-31G(2df,p)`.
- **`job.workflow.soc(...)`** is the one call that makes this a SOC run: it sets
  `runtype=soc`, forwards `soc_2e=1`, and applies `scal_rel=2` to the `[scf]`
  reference (the helper's default is already `2`).

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd spin-orbit-coupling/inputs
openqp h2o_soc.oqp
```

Python style:

```bash
cd spin-orbit-coupling/inputs
python h2o_soc.py
```

Both need OpenQP installed (`pip install openqp`) and compute the same SOC
matrix elements and eigenvalues.

## Reading the output

The workflow produces three things you will want:

- **The spin-free MRSF energies** — the singlet and triplet roots *before* SOC
  mixing. In Python these come back as `td_singlet_energies` and
  `td_triplet_energies` inside `mol.get_results()`.
- **The SOC matrix elements** between those states — the couplings that drive
  intersystem crossing and phosphorescence.
- **The SOC eigenvalues** — energies of the spin-adiabatic states obtained by
  diagonalizing the spin-free-plus-SOC Hamiltonian, reported in cm⁻¹. From
  Python, `mol.get_soc()` returns exactly these.

From the **Python** side the two calls to look at are:

```python
mol.get_soc()          # SOC eigenvalues, cm^-1 (spin-adiabatic states)
mol.get_results()      # dict: td_singlet_energies, td_triplet_energies, soc, ...
```

`mol.get_results()` is the JSON-friendly summary that mirrors the
`<project>.json` file. In the **log file** (`h2o_soc.log`) look for the converged
ROHF reference energy, the block of MRSF singlet and triplet roots, and the
printed SOC matrix / eigenvalues. A large SOC element between a given
singlet-triplet pair is the fingerprint of a fast spin-forbidden channel between
those two states.

## Manual

- Spin-orbit coupling workflow (the `runtype=soc` contract, `soc_2e`,
  `scal_rel`, and the reported quantities):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/spin-orbit-coupling/>
- MRSF-TDDFT workflow (the triplet-ROHF reference and `[tdhf] type=mrsf` roots
  that SOC couples):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/mrsf-tddft/>
- Running OpenQP from Python (the `job.theory.mrsf(...)` /
  `job.workflow.soc(...)` idiom):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
