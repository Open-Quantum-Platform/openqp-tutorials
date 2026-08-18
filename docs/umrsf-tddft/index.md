# UMRSF-TDDFT: spin-flip excited states from a UHF triplet reference

**Mixed-Reference Spin-Flip TDDFT (MRSF-TDDFT)** gives balanced, low-spin-contaminated
singlet *and* triplet excited states — including the near-degenerate,
diradical-character situations (conical intersections, bond breaking, polyenes)
where ordinary TDDFT fails. Standard MRSF spin-flips out of a **high-spin ROHF**
triplet reference. **UMRSF** ("U" for unrestricted) does the same physics on an
**unrestricted (UHF)** triplet reference instead. You reach for it when you want
the MRSF target states but prefer — or need — a UHF reference determinant rather
than ROHF orbitals. This tutorial computes single-point UMRSF-TDDFT excitation
energies for **butadiene (C₄H₆)** at BHHLYP / 6-31G\*.

## A little theory

MRSF-TDDFT starts from a **triplet (M_s = +1) reference** and applies a
**spin-flip** response operator: flipping one β-spin excitation down mixes the two
open-shell reference configurations, so the resulting response states are proper
**singlets and triplets** built from a common set of orbitals. Mixing the *two*
reference determinants (the "mixed reference") is what removes the spin
contamination that plagues plain spin-flip TDDFT and delivers a balanced
description of states that are about to cross.

The only thing UMRSF changes is the **reference wavefunction**: instead of a
restricted-open-shell (ROHF) triplet, it uses an **unrestricted (UHF) triplet**,
and the response is solved with the unrestricted spin-flip kernel. In OpenQP the
switch is two coordinated keywords — `[scf] type=uhf` (the "U") and
`[tdhf] type=umrsf` — with the reference forced to high spin
(`[scf] multiplicity=3`). Everything downstream (the BHHLYP functional, the basis,
the requested number of roots) is identical to ordinary MRSF. For the derivation
and the full method family, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/c4h6_umrsf_energy.oqp`](inputs/c4h6_umrsf_energy.oqp) —
butadiene in 6-31G\*, a UHF triplet reference, an unrestricted-MRSF response,
energy only. Annotated:

```text
umrsf(nstate=5)/bhhlyp/6-31g*    # model(roots)/functional/basis
geom="""
C   -1.901080641    0.114577198    0.000000000
C   -0.574900598   -0.402237631    0.000000000
C    1.901080641   -0.114577198    0.000000000
C    0.574900598    0.402237631    0.000000000
H   -2.755248674   -0.545115910    0.000000000
H   -2.071562576    1.183086840    0.000000000
H    2.755248674    0.545115910    0.000000000
H   -0.441348756   -1.477941036    0.000000000
H    2.071562576   -1.183086840    0.000000000
H    0.441348756    1.477941036    0.000000000
"""
```

What each part does:

- **`umrsf`** is the single token that makes this "U"MRSF. It carries all four of
  the legacy choices at once: an **unrestricted (UHF)** reference determinant — the
  "U" — in the **high-spin triplet** state that spin-flip must start from, the
  **unrestricted mixed-reference spin-flip** response built on it, and the
  **singlet** target manifold. Changing the one token to `mrsf` recovers ordinary
  MRSF-TDDFT on an ROHF reference; the geometry, functional, basis, and root count
  stay exactly as they are.
- **`nstate=5`** requests five target roots. It is a model option, so it lives in
  the route parentheses.
- **`bhhlyp` / `6-31g*`** are the standard MRSF functional/basis pairing. Drop the
  functional component (`umrsf-tdhf/6-31g*`) for UMRSF-TDHF.
- **No driver line means `energy()`** — single-point excitation energies. Add
  `grad(S1)` or `opt(S1)` when you want derivatives.

MRSF vs UMRSF, at a glance:

| | route | reference the model implies |
| --- | --- | --- |
| MRSF-TDDFT | `mrsf(nstate=5)/bhhlyp/6-31g*` | ROHF, triplet |
| UMRSF-TDDFT | `umrsf(nstate=5)/bhhlyp/6-31g*` | UHF, triplet |

To collect **triplet** target states instead of singlets, name a triplet in the
driver — `energy(T0)` — rather than setting a target multiplicity by hand.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/c4h6_umrsf_energy.py`](inputs/c4h6_umrsf_energy.py). The compact
`job.theory.mrsf(...)` helper always sets an **ROHF** reference and `type=mrsf`, so
it cannot express the UHF/`umrsf` variant. Instead the script writes the same
`[input]` / `[scf]` / `[tdhf]` keywords the deck uses directly, via
`job.update(...)` — the documented raw-section style for calculations the
model-specific helpers do not cover.

```python
from oqp.openqp import OpenQP

# Butadiene geometry (Angstrom), same atoms/order as the .oqp deck.
butadiene = """
6   -1.901080641    0.114577198    0.000000000
6   -0.574900598   -0.402237631    0.000000000
6    1.901080641   -0.114577198    0.000000000
6    0.574900598    0.402237631    0.000000000
1   -2.755248674   -0.545115910    0.000000000
1   -2.071562576    1.183086840    0.000000000
1    2.755248674    0.545115910    0.000000000
1   -0.441348756   -1.477941036    0.000000000
1    2.071562576   -1.183086840    0.000000000
1    0.441348756    1.477941036    0.000000000
"""

job = OpenQP("c4h6_umrsf", silent=1)
job.molecule(butadiene, charge=0)

# UMRSF-TDDFT = method=tdhf + UHF triplet reference + [tdhf] type=umrsf.
job.update({
    "input": {"method": "tdhf", "functional": "bhhlyp",
              "basis": "6-31g*", "runtype": "energy"},
    "scf":   {"type": "uhf", "multiplicity": 3},                 # UHF high-spin reference
    "tdhf":  {"type": "umrsf", "nstate": 5, "multiplicity": 1},  # 5 singlet targets
})

# Serial run for an ordinary (non-MPI) Python invocation.
job.control(usempi=False)

mol = job.run()
results = mol.get_results()

print("Reference (SCF) energy:", results["energy"])
print("UMRSF state energies (Hartree, relative to reference):")
for i, e in enumerate(results["td_energies"]):
    print(f"  state {i}: {e:+.6f}")
```

Each dictionary passed to `job.update(...)` maps one-to-one onto a `[section]` of
the input deck, so the three entries here are exactly the `[input]`, `[scf]`, and
`[tdhf]` blocks above. `job.control(usempi=False)` runs serially, which is what you
want for a plain `python …` invocation. To recover ordinary MRSF you would change
two values — `"scf": {"type": "rohf", …}` and `"tdhf": {"type": "mrsf", …}`.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd umrsf-tddft/inputs
openqp c4h6_umrsf_energy.oqp
```

Python style:

```bash
cd umrsf-tddft/inputs
python c4h6_umrsf_energy.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

A UMRSF-TDDFT single point produces one **reference (SCF) energy** plus a list of
**excited-state energies**, one per requested root:

- In the **log file** (`<project>.log`) look for the converged SCF energy of the
  UHF triplet reference, then the table of `nstate` MRSF roots with their
  excitation energies (and oscillator strengths).
- From **Python**, `mol.get_results()` returns a dict: `results["energy"]` is the
  reference energy and `results["td_energies"]` is the list of state energies the
  script iterates over. (`mol.get_scf_energy()` returns the same SCF reference
  energy.)
- With `[tdhf] multiplicity=1` the five roots are **singlet** target states; the
  gaps between successive `td_energies` are the singlet excitation energies of
  butadiene. Re-run with `multiplicity=3` to obtain the triplet manifold from the
  same reference.

## Manual

- MRSF-TDDFT / TDHF response (method background, the `[tdhf]` keyword contract, and
  the `mrsf` vs. `umrsf` reference choice):
  <https://open-quantum-platform.github.io/openqp-docs/>
- `[scf]` keyword reference (`type=uhf|rohf`, `multiplicity` — choosing the
  high-spin reference the spin flip starts from):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
