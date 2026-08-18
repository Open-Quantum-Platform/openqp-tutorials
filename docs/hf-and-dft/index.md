# Hartree-Fock and DFT energies and gradients

The self-consistent-field (SCF) energy is the foundation of nearly every
quantum-chemistry calculation. **Hartree-Fock (HF)** solves the electronic
Schrodinger equation in a mean-field approximation — each electron moves in the
averaged field of all the others — and **Kohn-Sham DFT** replaces HF's exact
exchange with an approximate exchange-correlation functional to fold in electron
correlation cheaply. This tutorial is the "hello world" of OpenQP: it computes a
single-point HF energy for water, then shows how to switch the reference
(RHF/ROHF/UHF), swap HF for DFT, and ask for a **nuclear gradient** (the forces
you need for geometry optimizations and dynamics). Everything downstream — MP2,
TDDFT, NAMD — starts from exactly the SCF step you build here.

## A little theory

An SCF calculation expands the molecular orbitals in your **basis set**
(here `6-31G*`), then iterates the Fock/Kohn-Sham equations until the orbitals
that build the potential match the orbitals the potential produces — the
"self-consistent" point. The **reference** determines how spin is handled:

- **RHF** (restricted) — closed-shell singlets; every spatial orbital holds a
  paired α/β electron. The right choice for water's ground state.
- **ROHF** (restricted open-shell) — open-shell systems kept as a spin
  eigenfunction; unpaired electrons share the same spatial orbitals as the paired
  ones.
- **UHF** (unrestricted) — α and β electrons get independent spatial orbitals;
  more variational freedom for radicals and bond-breaking, at the cost of spin
  contamination.

**DFT** uses the identical machinery; the only change is that `functional`
(e.g. `bhhlyp`) is non-empty, so the exchange-correlation potential replaces
pure HF exchange. A **gradient** run does one extra thing after the SCF
converges: it differentiates the energy with respect to the nuclear coordinates
to get the force on every atom. For the full derivation and the list of
functionals and references, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_rhf_energy.oqp`](inputs/h2o_rhf_energy.oqp) —
water in the 6-31G\* basis, a closed-shell RHF reference, single-point energy.
The whole calculation is four lines:

```text
rhf/6-31g*
guess(type=huckel)
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Key points:

- **`rhf/6-31g*`** is the *route*: model first, then the basis. `rhf` is pure
  Hartree-Fock with a restricted closed-shell reference — there is no separate
  "method" and "reference" to keep in sync. Use `uhf` or `rohf` for the other two
  references, and add `mult=3` for a triplet.
- To run **DFT** instead, put a functional between the model and the basis:
  `rks/bhhlyp/6-31g*`. That third component is exactly what turns the SCF into
  Kohn-Sham DFT; `uks` and `roks` are the unrestricted and restricted-open-shell
  Kohn-Sham models.
- **No driver line means `energy()`** — a single-point energy. To get forces, add
  a driver line:

  ```text
  grad
  ```

  (`grad(S0)` if you prefer to name the state explicitly; for a ground-state model
  they mean the same thing.)
- **`guess(type=huckel)`** is an *exact section call*: it sets `type` in the
  legacy `[guess]` section, seeding the SCF with cheap extended-Huckel orbitals.
  Anything the concise surface does not name has this escape hatch — for example
  `scf(conv=1e-10)` or `scf(stability=true)` for a UHF reference.
- **`geom`** takes either an inline coordinate block in triple quotes, as here, or
  a path: `geom="water.xyz"`. Element symbols and atomic numbers are both fine.

The four variants this tutorial discusses, as complete routes:

| Calculation | Route |
| --- | --- |
| RHF energy | `rhf/6-31g*` |
| ROHF triplet energy | `rohf/6-31g*` + `mult=3` |
| UHF triplet energy | `uhf/6-31g*` + `mult=3` + `scf(stability=true)` |
| BHHLYP DFT energy | `rks/bhhlyp/6-31g*` |
| RHF gradient | `rhf/6-31g*` + `grad` |

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_rhf_energy.py`](inputs/h2o_rhf_energy.py). `job.theory.hf(...)` sets
`method=hf` and the `[scf]` reference in one call; `job.theory.dft(...)` does the
same but fills in the `functional`; and `job.workflow.gradient(...)` replaces the
energy workflow with a gradient run. The first block reproduces the `.oqp`
exactly; the rest demonstrate the other references, DFT, and a gradient:

```python
from oqp.openqp import OpenQP

WATER = """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""

# --- RHF/6-31G* energy: the exact equivalent of h2o_rhf_energy.oqp ----------
job = OpenQP("h2o_rhf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.hf(reference="rhf", basis="6-31g*")   # closed-shell HF
mol = job.run()
print("RHF energy:", mol.get_scf_energy())        # -> -76.0107 Hartree

# --- ROHF triplet reference (open-shell, restricted) ------------------------
job = OpenQP("h2o_rohf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=3)
job.theory.hf(reference="rohf", basis="6-31g*")
print("ROHF energy:", job.run().get_scf_energy())

# --- UHF triplet reference (open-shell, unrestricted) -----------------------
job = OpenQP("h2o_uhf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=3)
job.theory.hf(reference="uhf", basis="6-31g*", stability=True)
print("UHF energy:", job.run().get_scf_energy())

# --- DFT (BHHLYP) energy: job.theory.dft instead of job.theory.hf -----------
job = OpenQP("h2o_dft_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", reference="rhf", basis="6-31g*")
print("DFT (BHHLYP) energy:", job.run().get_scf_energy())

# --- RHF gradient: swap the energy workflow for job.workflow.gradient -------
# state=0 is the SCF (ground-state) gradient; it maps to [properties] grad=0.
job = OpenQP("h2o_rhf_grad", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.hf(reference="rhf", basis="6-31g*")
job.workflow.gradient(state=0)
mol = job.run()
print("RHF energy:", mol.get_scf_energy())
print("RHF gradient (Hartree/Bohr):", mol.get_grad())
```

What each helper maps to:

- **`job.theory.hf(reference=...)`** → `[input] method=hf` + `[scf] type=...`
  (empty functional, so pure HF). `reference` is `rhf` / `rohf` / `uhf`.
- **`job.theory.dft(functional=..., reference=...)`** → same, but with the
  `functional` set, which switches the SCF to Kohn-Sham DFT.
- **`stability=True`** (UHF block) requests an SCF wavefunction-stability check —
  useful for open-shell/unrestricted references that can land in a saddle point.
- **`job.workflow.gradient(state=0)`** → `[properties] grad=0`, the ground-state
  SCF gradient. (Excited-state gradients use `state=1, 2, …`.) Without it, `run()`
  does the energy-only workflow.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd hf-and-dft/inputs
openqp h2o_rhf_energy.oqp
```

Python style (runs the RHF energy plus the ROHF/UHF/DFT/gradient demos):

```bash
cd hf-and-dft/inputs
python h2o_rhf_energy.py
```

Both need OpenQP installed (`pip install openqp`). The `.oqp` and the first
Python block produce the same RHF energy.

## Reading the output

The SCF converges to a single number — the **total electronic energy** — plus,
for a gradient run, a force per atom. For this water / 6-31G* example the RHF
energy is about **-76.0107 Hartree**.

- In the **log file** (`<project>.log`) look for the converged SCF/HF total
  energy, printed once the iterations reach the `conv` threshold. A gradient run
  additionally prints the Cartesian gradient (Hartree/Bohr) for each atom.
- From **Python**, `mol.get_scf_energy()` returns the converged SCF energy for
  any reference (RHF, ROHF, UHF) or the DFT total; `mol.get_grad()` returns the
  nuclear gradient array after a `job.workflow.gradient(...)` run.
- Switching the **reference** (RHF → ROHF → UHF) or **method** (HF → DFT) changes
  the energy but is read back through the same `get_scf_energy()` call — the
  triplet ROHF/UHF energies here are higher than the closed-shell RHF singlet, as
  expected for water.

## Manual

- `[scf]` keyword reference (`type`, `multiplicity` — choosing the reference and
  functional): <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- Running OpenQP from Python (the `job.theory.hf(...)` / `job.theory.dft(...)` /
  `job.workflow.gradient(...)` idioms):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
