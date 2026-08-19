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
choice of `c_ss` and `c_os`. In OpenQP, MP2 is a **post-SCF** workflow: it first
converges an HF reference, then adds the correlation energy and reports the total
`E(HF+MP2)`. An analytic **ground-state gradient** is available on an RHF
reference, so `mp2/6-31g` also supports `grad` and `opt`; UHF/ROHF MP2 is
energy-only. The reference can be closed-shell **RHF**,
unrestricted **UHF**, or restricted-open-shell **ROHF** (ROHF orbitals are
semicanonicalized before the correlation step so the energy denominators are well
defined). For the derivation and the parameterizations behind each preset, see the
[MP2 workflow page](https://open-quantum-platform.github.io/openqp-docs/workflows/mp2/).

## Input-file style

The runnable deck is [`inputs/h2o_mp2_6-31g.oqp`](inputs/h2o_mp2_6-31g.oqp) —
water in the 6-31G basis, a UHF reference, conventional MP2. Annotated:

```text
mp2(reference=uhf)/6-31g      # model(reference)/basis
scf(conv=1e-10)               # tight SCF so the correlation energy is clean
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Key points:

- **`mp2/6-31g`** is `model/basis` — MP2 takes **no functional component**, and the
  route rejects one, because MP2 needs a Hartree-Fock reference rather than a
  Kohn-Sham one. That whole class of mistake is unavailable rather than merely
  discouraged.
- **`reference=uhf`** picks the HF reference: `rhf`, `uhf`, or `rohf`. Water is
  closed-shell, so `rhf` would be the natural choice and gives the identical energy
  here; this example ships with `uhf` to exercise the unrestricted path. For an
  open-shell system add `mult=3` (or whatever applies) as a top-level option.
- **`variant`** selects the spin-scaling preset; conventional MP2 is the default,
  so it is not written here. To use spin-component scaling, add the one route
  option:

  ```text
  mp2(reference=rhf,variant=scs-mp2)/6-31g    # Grimme SCS-MP2: c_ss = 1/3, c_os = 1.2
  ```

  ```text
  mp2(reference=rhf,variant=sos-mp2)/6-31g    # scaled-opposite-spin: c_ss = 0, c_os = 1.3
  ```

  Other accepted values include `os-mp2`, `ss-mp2`, `scs-mi-mp2`, and `custom`. For
  a literature parameterization not in the table, use `custom` with explicit scales:

  ```text
  mp2(reference=rhf,variant=custom,same_spin_scale=0.50,opposite_spin_scale=1.10)/6-31g
  ```

- **No driver line means `energy()`**. MP2 also has an analytic **ground-state**
  gradient, but only on an **RHF** reference: `mp2/6-31g` + `grad` works, while
  asking for a gradient on the `uhf` reference used here is rejected rather than
  silently downgraded. Excited-state MP2 gradients do not exist.
- **`scf(conv=1e-10)`** tightens the reference convergence — correlation energies
  are sensitive to it. `scf(...)` is an exact call into the legacy `[scf]`
  section, the escape hatch for anything the concise surface does not name; the
  initial guess, iteration cap, and so on stay at their defaults and are not
  written.

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
openqp h2o_mp2_6-31g.oqp
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
