# MP2 and spin-scaled MP2: adding correlation on top of Hartree-Fock

Hartree-Fock (HF) gives you a mean-field wavefunction: every electron feels the
*average* field of all the others, but their instantaneous **correlation** — the
way electrons dodge each other — is missing. **Second-order Moller-Plesset
perturbation theory (MP2)** is the cheapest, most widely used way to put that
correlation back. You run it when HF energies aren't accurate enough (reaction
energies, conformers, weak interactions) but a full coupled-cluster treatment is
too expensive. This tutorial computes a standalone MP2 energy for water and shows
how to switch on the **spin-component-scaled** variants (SCS-MP2, SOS-MP2) and how
the choice of HF reference (RHF/UHF/ROHF) enters.

## A little theory

MP2 treats electron correlation as a perturbation on the converged HF reference.
Its correlation energy splits cleanly into two physically distinct pieces — the
**same-spin** part (two electrons of the same spin, `E_aa + E_bb`) and the
**opposite-spin** part (`E_ab`):

```text
E(MP2) = c_ss * (E_aa + E_bb) + c_os * E_ab
```

Conventional MP2 sets both scale factors to 1. But MP2 systematically
*over*-counts same-spin correlation and *under*-counts opposite-spin correlation,
so re-weighting the two pieces improves accuracy for almost no extra cost. That is
the idea behind the **spin-component-scaled** family:

- **SCS-MP2** (Grimme): `c_ss = 1/3`, `c_os = 1.2` — a balanced all-round
  improvement.
- **SOS-MP2** (Head-Gordon): `c_ss = 0`, `c_os = 1.3` — drops same-spin work
  entirely (cheaper) and is well suited to Laplace/RI accelerations.

OpenQP computes the two spin components separately, so every preset is just a
choice of `c_ss` and `c_os`. In OpenQP, MP2 is a **post-SCF, energy-only**
workflow: it first converges an HF reference, then adds the correlation energy and
reports the total `E(HF+MP2)`. The reference can be closed-shell **RHF**,
unrestricted **UHF**, or restricted-open-shell **ROHF** (ROHF orbitals are
semicanonicalized before the correlation step so the energy denominators are well
defined). For the derivation and the parameterizations behind each preset, see the
[MP2 workflow page](https://open-quantum-platform.github.io/openqp-docs/workflows/mp2/).

## Input-file style

The runnable deck is [`inputs/h2o_mp2_6-31g.inp`](inputs/h2o_mp2_6-31g.inp) —
water in the 6-31G basis, a UHF reference, conventional MP2. Annotated:

```ini
[input]
system=
 8   0.000000000   0.000000000  -0.041061554   # O   (Angstrom)
 1  -0.533194329   0.533194329  -0.614469223   # H
 1   0.533194329  -0.533194329  -0.614469223   # H

charge=0
method=mp2              # select the MP2 post-SCF workflow
basis=6-31g
runtype=energy          # MP2 is energy-only; other runtypes are rejected
functional=             # MUST be empty — MP2 needs an HF reference, not KS-DFT

[guess]
type=huckel             # extended-Huckel initial orbital guess
save_mol=false

[scf]
type=uhf                # the HF reference: rhf | uhf | rohf
multiplicity=1          # singlet ground state
maxit=50                # max SCF iterations
conv=1.0e-10            # tight SCF convergence (correlation is sensitive)
save_molden=false

[mp2]
variant=mp2             # conventional (unscaled) MP2; c_ss = c_os = 1.0
```

Key points:

- **`method=mp2`** in `[input]` is what turns on the correlation step; **`functional`
  must be empty** (MP2 on a Kohn-Sham reference is rejected) and **`runtype` must be
  `energy`** (no MP2 gradients/Hessians yet).
- The **HF reference is chosen in `[scf]`** via `type`. Water is closed-shell, so
  `type=rhf` would be the natural choice and gives the identical energy here; this
  example ships with `uhf` to exercise the unrestricted path. Use `rohf` (or `uhf`)
  with the appropriate `multiplicity` for open-shell systems.
- The **`[mp2]` section is optional**. Omit it entirely for conventional MP2, or set
  `variant` to pick a spin-scaled preset. For spin-component scaling change one line:

  ```ini
  [mp2]
  variant=scs-mp2         # Grimme SCS-MP2: c_ss = 1/3, c_os = 1.2
  ```

  ```ini
  [mp2]
  variant=sos-mp2         # scaled-opposite-spin: c_ss = 0, c_os = 1.3
  ```

  Other accepted values include `os-mp2`, `ss-mp2`, `scs-mi-mp2`, and `custom`. For
  a literature parameterization not in the table, use `custom` with explicit scales:

  ```ini
  [mp2]
  variant=custom
  same_spin_scale=0.50
  opposite_spin_scale=1.10
  ```

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_mp2_6-31g.py`](inputs/h2o_mp2_6-31g.py). `job.theory.mp2(...)` sets
the `[input] method=mp2`, the `[scf]` reference, and the `[mp2]` section in one
call; extra keywords (here `conv`) are forwarded to `[scf]`.

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_mp2", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# reference -> [scf] type ; variant -> [mp2] variant ; conv -> [scf] conv
job.theory.mp2(basis="6-31g", reference="uhf", variant="mp2", conv=1.0e-10)

mol = job.run()
results = mol.get_results()
print("HF reference energy:", mol.get_scf_energy())
print("MP2 total energy:   ", results["energy"])
```

To run spin-scaled MP2, change only the `variant` argument:

```python
job.theory.mp2(basis="6-31g", reference="uhf", variant="scs-mp2")   # SCS-MP2
job.theory.mp2(basis="6-31g", reference="uhf", variant="sos-mp2")   # SOS-MP2
```

Custom scales go through the same helper (the helper defaults `variant` to
`custom` when you pass scale factors):

```python
job.theory.mp2(
    basis="6-31g",
    reference="uhf",
    variant="custom",
    same_spin_scale=0.50,
    opposite_spin_scale=1.10,
)
```

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd mp2/inputs
openqp h2o_mp2_6-31g.inp
```

Python style:

```bash
cd mp2/inputs
python h2o_mp2_6-31g.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

MP2 reports an **HF reference energy**, the **MP2 correlation energy** it adds, and
their sum, the **MP2 total energy** — the number you almost always want. For this
water / 6-31G example the validated run gives:

| Quantity | Value (Ha) |
| --- | --- |
| `E(MP2, correlation)` | `-0.1278307451` |
| `E(MP2, total)` | `-76.1121207760` |

- In the **log file** (`<project>.log`) look for the HF energy, the printed
  correlation energy, and the final total.
- From **Python**, `mol.get_scf_energy()` is the HF reference and
  `mol.get_results()["energy"]` is the MP2 total (this matches the `energy` field
  written to `<project>.json`).
- A **spin-scaled** run reuses the same HF reference but reweights the two
  correlation components, so its total differs from conventional MP2 by the change
  in `c_ss (E_aa+E_bb) + c_os E_ab` — the HF energy is unchanged.

## References / manual

- MP2 workflow (full input contract, spin-scaling table, implementation notes):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/mp2/>
- `[mp2]` keyword reference (`variant`, `same_spin_scale`, `opposite_spin_scale`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/mp2/>
- `[scf]` keyword reference (`type`, `multiplicity`, `conv` — choosing the HF
  reference): <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- Running OpenQP from Python (the `job.theory.mp2(...)` idiom):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
