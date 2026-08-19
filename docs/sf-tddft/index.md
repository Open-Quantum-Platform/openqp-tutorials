# Spin-flip TDDFT: excited states and diradicals from a triplet reference

Ordinary linear-response TDDFT excites electrons out of a **closed-shell** ground
state, which works well for bright, single-excitation valence states but breaks
down exactly where chemistry gets interesting: **conical intersections**, bond
dissociation, and **diradicals**, where two configurations become degenerate and
a single-reference ground state is qualitatively wrong. **Spin-flip TDDFT
(SF-TDDFT)** sidesteps that failure. Instead of exciting from the singlet, it
starts from a well-behaved **high-spin (triplet) reference** and reaches the
singlet manifold — including the multiconfigurational ground state — through
**spin-flipping** excitations. You reach for SF-TDDFT when you need a balanced
description of near-degenerate states that plain TDDFT cannot give. This tutorial
computes the SF-TDDFT roots of water and the gradient of the lowest root.

## A little theory

A spin-flip excitation takes the triplet M_s = +1 reference (two unpaired
same-spin electrons) and flips **one** spin down (α → β). The resulting
determinants are M_s = 0, so they span the **singlet** states — and, because you
generate them all on the same footing from a single high-spin reference, states
that are single-reference and states that are multireference are treated with the
**same** balanced accuracy. That is the whole trick: the problematic
near-degeneracy in the target singlet manifold is handled by a reference that is
itself far from degenerate.

Two consequences show up directly in the deck:

- The reference must be **high-spin ROHF** (`[scf] multiplicity=3 type=rohf`), not
  the usual closed-shell RHF. It is a mathematical starting point, not the state
  you care about.
- The SF **roots are re-indexed**: root 1 is the **lowest** spin-flip state (here
  the singlet ground state, S₀), root 2 the first excited singlet, and so on. So
  "the ground-state gradient" is the gradient of SF **root 1**, and higher roots
  are excited states.

A **half-and-half** functional such as **BHHLYP** (50 % exact exchange) is the
standard recommendation for spin-flip: the large fraction of Hartree-Fock
exchange keeps the spin-flip block well-behaved and reduces spin contamination.
For the full derivation and the related mixed-reference variant (MRSF-TDDFT, which
removes the residual spin contamination of plain SF), see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_sf-tddft.oqp`](inputs/h2o_sf-tddft.oqp) — water
in the 6-31G\* basis, BHHLYP, energy **and** gradient of an SF root. Annotated:

```text
sf(nstate=3)/bhhlyp/6-31g* grad(root=3)    # model(roots)/functional/basis + gradient of SF root 3
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Key points:

- **`sf(nstate=3)/bhhlyp/6-31g*`** is the whole electronic model. `sf` carries the
  **high-spin triplet ROHF reference** that spin-flip requires — that reference is
  non-negotiable for SF, so the format supplies it rather than asking you to.
  `nstate=3` solves the three lowest spin-flip roots. `bhhlyp` is the recommended
  half-and-half functional; `sf-tdhf/6-31g*` drops the functional for SF-TDHF.
- **`grad(root=3)`** is the driver. Note the `root=` spelling: unlike MRSF, the
  spin-flip roots are **not spin-adapted before diagonalization**, so their
  physical labels are not known in advance and there is no honest `S1` to write.
  `root=1` is the lowest (multiconfigurational ground) state, so `root=3` is the
  gradient of the second excited SF root. Use `energy()` for energies only.

This is the one place the format asks for an implementation index instead of a
physical label — and it does so because the physics, not the interface, is what
makes the label unavailable.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_sf-tddft.py`](inputs/h2o_sf-tddft.py). `job.theory.sf_tddft(...)`
sets `[tdhf] type=sf` and defaults to the high-spin reference
(`reference="rohf"`, `multiplicity=3`) that SF requires, so you do not spell out
the `[scf]` block; `job.workflow.gradient(state=3)` selects `runtype=grad` and
`[properties] grad=3`.

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_sf_tddft", silent=1)

# Water geometry (Angstrom), same as the input deck.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
)

# Spin-flip TDDFT: the helper sets [tdhf] type=sf and defaults to the
# high-spin reference (reference="rohf", multiplicity=3) that SF requires.
job.theory.sf_tddft(functional="bhhlyp", basis="6-31g*", nstate=3)

# Gradient of SF root 3. For SF/MRSF, state 1 is the lowest root (the
# multiconfigurational ground state); this maps to [properties] grad=3.
job.workflow.gradient(state=3)

mol = job.run()

results = mol.get_results()
print("Reference (ROHF) energy:", results["energy"])
print("SF-TDDFT root energies:", results["td_energies"])
print("Gradient (root 3):", mol.get_grad())
```

The two calls carry all the physics: `job.theory.sf_tddft(...)` fixes the
functional, basis, spin-flip type, and root count, and `job.workflow.gradient(state=3)`
fixes the run type and which root is differentiated. Everything else (the triplet
ROHF reference, the Hückel guess) is filled in by the SF defaults.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd sf-tddft/inputs
openqp h2o_sf-tddft.oqp
```

Python style:

```bash
cd sf-tddft/inputs
python h2o_sf-tddft.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

The calculation reports three things:

- the **reference (ROHF) energy** — the high-spin triplet the SF states are built
  from,
- the **SF-TDDFT root energies** — the spin-flip manifold, ordered lowest-first
  (root 1 = S₀, root 2 = S₁, …), and
- the **gradient** of the requested root (root 3 here).

Where to find them:

- In the **log file** (`<project>.log`) look for the ROHF energy, the table of
  spin-flip excitation energies / total root energies, and the gradient block.
- From **Python**, `mol.get_results()["energy"]` is the reference (ROHF) energy,
  `mol.get_results()["td_energies"]` is the list of SF root energies, and
  `mol.get_grad()` returns the gradient of the root selected by
  `job.workflow.gradient(state=...)`.
- Because the SF roots are re-indexed, the **first** entry of `td_energies` is the
  ground state; higher entries are the excited singlets. Match `grad`'s index to
  the root you actually want a force for (`grad=1` for the ground-state gradient).

## Manual

- Spin-flip / MRSF-TDDFT and the `[tdhf]` response engine (full input contract,
  `type=sf` vs `mrsf`, root indexing):
  <https://open-quantum-platform.github.io/openqp-docs/>
- `[scf]` keyword reference (`type`, `multiplicity` — choosing the high-spin
  reference SF requires):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
