# TDDFT and TDHF excited states: vertical excitation energies

Ground-state Hartree-Fock (HF) and DFT tell you where a molecule sits in its
lowest electronic state, but photochemistry, UV-Vis spectra, and fluorescence all
live in the **excited** states. **Linear-response time-dependent** theory is the
workhorse for getting them: **TDDFT** (time-dependent DFT) and its wavefunction
cousin **TDHF** (time-dependent Hartree-Fock, a.k.a. the random-phase
approximation) give you vertical excitation energies — the energy to promote the
molecule from its ground state to each excited state *at the fixed ground-state
geometry*. You run these when you want an absorption spectrum, the character of an
excited state, or a cheap first look at photophysics before committing to
gradients or dynamics. This tutorial computes the three lowest singlet excitation
energies of water and shows how one keyword switches between TDDFT and TDHF.

## A little theory

Both methods answer the same question — "what are the excited states as a response
of the ground state to a perturbation?" — by solving the same
**Casida / random-phase-approximation** eigenvalue problem built from the occupied
and virtual orbitals of a converged reference. The eigenvalues are the excitation
energies; the eigenvectors are the excited-state wavefunctions (which orbital
transitions dominate). The **only** difference is the reference and the response
kernel:

- **TDHF** uses a Hartree-Fock reference and the exact-exchange kernel — pure
  wavefunction linear response, no functional.
- **TDDFT** uses a Kohn-Sham (DFT) reference and adds the exchange-correlation
  kernel of the chosen functional on top. In practice TDDFT is more accurate and
  much more widely used because the functional recovers correlation the HF
  reference lacks.

In OpenQP a single driver, `method=tdhf`, runs *both*: leave the `functional`
empty and you get TDHF; set a functional and the same driver becomes TDDFT. This
tutorial uses the hybrid **B3LYP5** functional and requests three roots. For the
derivation and the full keyword contract, see the
[TDDFT/TDHF workflow page](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is
[`inputs/h2o_tddft_energy.inp`](inputs/h2o_tddft_energy.inp) — water in the 6-31G\*
basis, a closed-shell RHF reference, B3LYP5 TDDFT, three excited-state roots.
Annotated:

```ini
[input]
system=
   8   0.000000000   0.000000000  -0.041061554   # O   (Angstrom)
   1  -0.533194329   0.533194329  -0.614469223   # H
   1   0.533194329  -0.533194329  -0.614469223   # H
charge=0
runtype=energy          # vertical excitation energies at a fixed geometry
basis=6-31g*
functional=b3lyp5       # non-empty -> TDDFT; leave empty for plain TDHF
method=tdhf             # the linear-response driver (drives BOTH TDDFT and TDHF)

[guess]
type=huckel             # extended-Huckel initial orbital guess

[scf]
multiplicity=1          # closed-shell singlet ground state
type=rhf                # restricted HF/KS reference

[tdhf]
nstate=3                # number of excited-state roots to solve
```

Key points:

- **`method=tdhf`** in `[input]` selects the linear-response driver. It is the
  common engine for *both* models — the presence or absence of **`functional`** is
  what distinguishes them. With `functional=b3lyp5` you get **TDDFT**; delete the
  `functional` line (or leave it empty) and the same deck becomes plain **TDHF**.
- **`runtype=energy`** requests excitation energies at the input geometry (vertical
  excitations). No excited-state gradient is computed here.
- The **reference is a closed-shell singlet**, so `[scf] type=rhf` and
  `multiplicity=1`. The response is built on top of this converged reference.
- **`[tdhf] nstate=3`** asks for the three lowest roots. Raise it to resolve more
  of the spectrum.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_tddft_energy.py`](inputs/h2o_tddft_energy.py). `job.theory.tddft(...)`
sets `method=tdhf` *with* a functional (making it TDDFT) on an RHF reference, and
`nstate` is the number of roots.

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_tddft", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# TDDFT: RHF reference + B3LYP5 functional, 3 excited-state roots.
job.theory.tddft(functional="b3lyp5", basis="6-31g*", nstate=3)

mol = job.run()
results = mol.get_results()

print("Ground-state (RHF) energy [Ha]:", results["energy"])
print("Excitation energies [Ha]     :", results["td_energies"])
```

For **plain TDHF**, swap the theory line for the no-functional variant (Hartree-Fock
linear response):

```python
job.theory.tdhf(reference="rhf", basis="6-31g*", nstate=3)
```

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd tddft-and-tdhf/inputs
openqp h2o_tddft_energy.inp
```

Python style:

```bash
cd tddft-and-tdhf/inputs
python h2o_tddft_energy.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

A TDDFT/TDHF energy run reports the **converged ground-state (SCF) energy** and one
**vertical excitation energy per requested root** — the numbers you almost always
want for a spectrum. The excitation energies come back in Hartree (multiply by
27.2114 for eV).

- In the **log file** (`<project>.log`) look for the ground-state energy followed by
  the table of excited states — for each root its excitation energy and dominant
  orbital transition.
- From **Python**, `results["energy"]` (equivalently `mol.get_scf_energy()`) is the
  ground-state RHF/KS energy, and `results["td_energies"]` is the list of the three
  excitation energies — exactly the two quantities the script prints.
- **Switching to TDHF** (empty functional) reuses the same machinery but with the
  HF reference and exact-exchange kernel, so both the ground-state energy and the
  excitation energies change; the number of roots and the way you read them are
  identical.

## References / manual

- OpenQP manual (workflows, keyword reference, and the `job.theory.tddft(...)` /
  `job.theory.tdhf(...)` idioms):
  <https://open-quantum-platform.github.io/openqp-docs/>
- `[tdhf]` keyword reference (`nstate`, and the shared TDDFT/TDHF driver):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/tdhf/>
- `[scf]` keyword reference (`type`, `multiplicity` — choosing the reference):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- Running OpenQP from Python:
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
