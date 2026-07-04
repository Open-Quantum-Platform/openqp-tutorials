# NMR shielding (GIAO / CGO): nuclear magnetic shielding tensors from an SCF wavefunction

Nuclear magnetic resonance (NMR) chemical shifts are one of chemistry's most
information-rich observables — they report on the electronic environment of every
magnetic nucleus in a molecule. The quantity a quantum-chemistry code computes is
the **nuclear magnetic shielding tensor** σ, which links the external magnetic
field to the field actually felt at a nucleus; the experimental chemical shift is
just the shielding measured relative to a reference compound. You run this
calculation to assign spectra, predict shifts for a proposed structure, or
disentangle overlapping resonances. This tutorial computes the shielding tensors
of water at the RHF/STO-3G level using **gauge-including atomic orbitals (GIAO)**,
and shows the same run through the OpenQP Python API.

## A little theory

A magnetic shielding calculation is a **second-derivative (response) property**:
σ is the mixed second derivative of the energy with respect to the external
magnetic field **B** and each nuclear magnetic moment **m**. The awkward part is
that the vector potential of a uniform field, **A = ½ B × (r − R₀)**, depends on
an arbitrary **gauge origin R₀**. With a finite (incomplete) basis set, a naive
calculation gives shieldings that *change* when you move R₀ — an unphysical
artifact.

Two schemes cure this:

- **GIAO** — *gauge-including* (a.k.a. *London*) *atomic orbitals* attach a
  field-dependent complex phase factor to every basis function, so each orbital
  carries its own local gauge origin. The result is **gauge-origin independent**
  and converges quickly with basis size. This is the modern default.
- **CGO** — a *common gauge origin* places a single R₀ for the whole molecule (the
  classical approach). It is simpler but gauge-dependent at finite basis; the
  answer depends on where you put the origin.

OpenQP selects between them with a single keyword, `nmr_gauge` (`giao` or `cgo`),
after an ordinary SCF wavefunction has converged. Shielding is a **post-SCF
property job**: converge the RHF (or DFT) reference, then solve the coupled-
perturbed response equations for the magnetic perturbation and assemble σ. For the
theory and the response formalism, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_giao-nmr.inp`](inputs/h2o_giao-nmr.inp) — water
in the minimal STO-3G basis, a closed-shell RHF reference, GIAO shielding.
Annotated:

```ini
[input]
system=
   8   0.000000000   0.000000000  -0.041061554   # O   (Angstrom)
   1  -0.533194329   0.533194329  -0.614469223   # H
   1   0.533194329  -0.533194329  -0.614469223   # H
charge=0
runtype=energy          # shielding is a post-SCF property on a single-point wavefunction
basis=sto-3g
method=hf               # Hartree-Fock reference (not KS-DFT here)

[guess]
type=huckel             # extended-Huckel initial orbital guess

[scf]
multiplicity=1          # closed-shell singlet
type=rhf                # restricted Hartree-Fock reference

[properties]
scf_prop=nmr            # request the NMR shielding property after the SCF
nmr_gauge=giao          # gauge treatment: giao (gauge-including) or cgo (common origin)
```

Key points:

- **`[input] method=hf`** with **`runtype=energy`** converges an ordinary
  single-point wavefunction; the shielding is layered on top afterwards, so no
  special run type is needed.
- The **`[scf]` section** fixes the reference: `type=rhf` and `multiplicity=1` for
  closed-shell water. The `[guess] type=huckel` line just seeds the SCF with an
  extended-Hückel guess.
- The **`[properties]` section is what turns the shielding on.** `scf_prop=nmr`
  asks OpenQP to evaluate the NMR shielding tensor as a post-SCF property, and
  `nmr_gauge` picks the gauge scheme. To switch from London orbitals to a common
  gauge origin, change one line:

  ```ini
  [properties]
  scf_prop=nmr
  nmr_gauge=cgo         # common gauge origin instead of GIAO
  ```

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_giao-nmr.py`](inputs/h2o_giao-nmr.py). `job.theory.hf(...)` sets the
HF reference and basis, and `job.workflow.nmr(gauge="giao")` maps onto
`[properties] scf_prop=nmr, nmr_gauge=giao`:

```python
from oqp.openqp import OpenQP

# Same water geometry as the input deck (Angstrom).
water = """
8   0.000000000   0.000000000  -0.041061554
1  -0.533194329   0.533194329  -0.614469223
1   0.533194329  -0.533194329  -0.614469223
"""

job = OpenQP("h2o_giao_nmr", silent=1)

# Molecular identity: closed-shell singlet, neutral.
job.molecule(water, charge=0, multiplicity=1)

# Quantum theory: Hartree-Fock reference SCF in the STO-3G basis (RHF by default).
job.theory.hf(basis="sto-3g")

# Workflow: NMR shielding with the gauge-including-atomic-orbital (GIAO) gauge.
# This maps onto [properties] scf_prop=nmr, nmr_gauge=giao.
job.workflow.nmr(gauge="giao")

mol = job.run()

results = mol.get_results()
print("SCF energy (Hartree):", results["energy"])
# Per-atom isotropic shielding is the last entry of each nmr_shielding row (ppm).
for atom, row in zip(results["atoms"], results["nmr_shielding"]):
    print(f"atom Z={atom:>2}  isotropic shielding = {row[-1]:10.4f} ppm")
```

To run the common-gauge-origin variant instead, change only the workflow call:

```python
job.workflow.nmr(gauge="cgo")   # common gauge origin instead of GIAO
```

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd nmr-shielding/inputs
openqp h2o_giao-nmr.inp
```

Python style:

```bash
cd nmr-shielding/inputs
python h2o_giao-nmr.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

The job first converges the SCF, then reports one **nuclear magnetic shielding
tensor per atom** (in ppm). The number chemists compare to experiment is the
**isotropic shielding**, the average of the tensor's three diagonal elements
(⅓·tr σ).

- From **Python**, `mol.get_results()` returns a dictionary. The script uses two
  fields: `results["energy"]` is the converged SCF energy (Hartree), and
  `results["nmr_shielding"]` is the per-atom shielding data — one row per atom,
  aligned with `results["atoms"]` (the nuclear charges Z). The script prints the
  **last entry of each row**, `row[-1]`, which is that atom's isotropic shielding
  in ppm.
- In the **log file** (`<project>.log`) look for the SCF total energy and the
  printed shielding tensor for each nucleus.
- Interpreting the numbers: a **larger** isotropic shielding means the nucleus is
  more shielded (electron density screening it from the field) and therefore
  appears **upfield** (smaller chemical shift). For water you will see the oxygen
  strongly shielded and the two symmetry-equivalent hydrogens with equal, much
  smaller shieldings. Absolute STO-3G values are far from quantitative — the point
  of this deck is the *workflow*; production shieldings need a larger basis (and
  usually a correlated or DFT reference).

## Manual

- NMR shielding / properties — the `scf_prop` and `nmr_gauge` keywords and the
  GIAO vs. CGO gauge treatments:
  <https://open-quantum-platform.github.io/openqp-docs/>
- Running OpenQP from Python (the `job.theory.hf(...)` / `job.workflow.nmr(...)`
  idiom): <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
