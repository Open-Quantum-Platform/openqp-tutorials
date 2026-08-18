# Hessians, frequencies, IR and Raman: the second derivatives of the energy

The energy gradient tells you the **forces** on the nuclei; the **Hessian** — the
matrix of second derivatives of the energy with respect to nuclear positions —
tells you the *curvature* of the potential energy surface. That single object
unlocks a lot: diagonalising the mass-weighted Hessian gives the **harmonic
vibrational frequencies** and **normal modes**, its eigenvalue signs confirm
whether a stationary point is a true minimum (all positive) or a transition state
(one negative), and combined with dipole and polarizability derivatives it yields
the **IR** and **Raman** intensities you compare against an experimental spectrum.
You run a Hessian to characterise an optimized geometry, to get zero-point energy
and thermochemistry, or to simulate a vibrational spectrum. This tutorial computes
an analytical HF/DFT Hessian for water and reads off the frequencies and spectra.

## A little theory

At a geometry **R** the energy expands as
`E(R + dR) = E(R) + g·dR + ½ dRᵀ H dR + …`, where `g` is the gradient and
`H = ∂²E/∂R∂R` is the Hessian. To get physical vibrations you **mass-weight** it,
`H̃_ij = H_ij / √(mᵢmⱼ)`, and diagonalise. Six eigenvalues (five for a linear
molecule) come out near zero — they are overall translation and rotation — and the
rest are `ω² = λ`, i.e. the eigenvalue is the square of a **harmonic frequency**
and its eigenvector is the corresponding **normal mode**. Each mode's **IR
intensity** is set by how much the dipole moment changes along it (`∂μ/∂Q`), and
its **Raman activity** by how much the polarizability changes (`∂α/∂Q`); OpenQP
evaluates both alongside the frequencies. Because the frequencies feed directly
into the vibrational partition function, the same run also prints **zero-point
energy** and finite-temperature **thermochemistry**. OpenQP can build `H` either
**analytically** (solving the coupled-perturbed equations — fast and accurate) or
by **numerical** finite differences of the gradient. For the details of both paths
see the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_hess.oqp`](inputs/h2o_hess.oqp) — water at a
BHHLYP/6-31G\* Kohn-Sham level, taking an analytical ground-state Hessian.
Annotated:

```text
rks/bhhlyp/6-31g*                     # restricted Kohn-Sham, BHHLYP, 6-31G*
hess(S0,type=analytical,clean=true)   # Hessian of the ground state
geom="""
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
"""
```

Key points:

- **`hess(...)` is the driver** that selects the Hessian workflow. It runs the SCF
  first, then builds the second-derivative matrix and diagonalises it for
  frequencies, normal modes, thermochemistry, and IR/Raman intensities — all
  written to the log.
- **`S0`** names the state physically: the ground state. Excited-state Hessians use
  the same spelling — `hess(S1)` on an `mrsf` route, for example.
- **`rks/bhhlyp/6-31g*`** is a Kohn-Sham DFT reference. Drop the functional
  component (`rhf/6-31g*`) for a pure Hartree-Fock Hessian; the `hess(...)`
  mechanics are identical either way. Use `uks`/`roks` with a `mult=` line for
  open-shell systems.
- **`type=analytical`** uses the coupled-perturbed (CPHF/CPKS) equations — the
  fast, accurate choice; switch to `type=numerical` to build the Hessian by finite
  differences of the gradient (useful when analytic second derivatives are not
  available for a method). `clean=true` removes the scratch Hessian files
  afterwards. Other `hess(...)` options include `dx` (displacement size for the
  numerical path), `nproc`, `restart`, and `temperature`.

> **Geometry note.** A Hessian is only meaningful at a **stationary point**. This
> water geometry is already optimized for BHHLYP/6-31G\*; if you compute a Hessian
> at an un-optimized geometry the residual forces contaminate the low frequencies.
> Optimize first (`opt(S0)`), then run the Hessian on that geometry.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_hess.py`](inputs/h2o_hess.py). `job.theory.dft(...)` sets the
KS-DFT reference (`method=hf` + `functional=bhhlyp` + `[scf] type=rhf`), and
`job.workflow.hessian(...)` selects `runtype=hess` and fills the `[hess]` section:

```python
from oqp.openqp import OpenQP

# silent=1 keeps the console quiet; the .log still records everything.
job = OpenQP("h2o_hess", silent=1)

# The built-in "water" geometry matches the deck; you can also pass an inline
# XYZ block or an .xyz path here.
job.molecule(geometry="water", charge=0, multiplicity=1)

# BHHLYP/6-31G* Kohn-Sham: method=hf + functional=bhhlyp + [scf] type=rhf.
job.theory.dft(functional="bhhlyp", basis="6-31g*")

# The analytical ground-state Hessian: runtype=hess + [hess] section.
# state=0 is the ground state; clean=True removes scratch Hessian files.
job.workflow.hessian(type="analytical", state=0, clean=True)

mol = job.run()

print("SCF energy (Hartree):", mol.get_scf_energy())

# Frequencies, normal modes, thermochemistry, and IR/Raman intensities are in
# the .log. The raw Cartesian Hessian is available here for post-processing.
hessian = mol.get_hess()
print("Hessian shape:", None if hessian is None else hessian.shape)

# get_results() is the JSON-friendly summary (atoms, coords, energy, hess, ...).
print(mol.get_results().keys())
```

For a **pure Hartree-Fock** Hessian instead of DFT, drop the functional and use the
HF theory helper (the `[hess]` call is unchanged):

```python
job.theory.hf(basis="6-31g*")                       # no functional -> HF reference
job.workflow.hessian(type="analytical", state=0)
```

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd vibrational-analysis/inputs
openqp h2o_hess.oqp
```

Python style:

```bash
cd vibrational-analysis/inputs
python h2o_hess.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

The Hessian run computes one SCF energy and then a whole table of vibrational
data. The main things to look at:

- **In the log file** (`<project>.log`): the **harmonic frequencies** (cm⁻¹) with
  their **IR intensities** and **Raman activities**, the **normal-mode**
  displacement vectors, and the **thermochemistry** block (zero-point energy,
  enthalpy, entropy, free energy). Water is bent and non-linear, so you get
  **3N − 6 = 3** real vibrational modes — the symmetric stretch, the asymmetric
  stretch, and the bend — plus six near-zero translation/rotation modes you can
  ignore. All three real frequencies should be **positive**, confirming this
  geometry is a genuine minimum; a single **imaginary** (reported negative)
  frequency would flag a transition state.
- **From Python**:
  - `mol.get_scf_energy()` — the converged SCF (here KS-DFT) energy in Hartree.
  - `mol.get_hess()` — the raw Cartesian Hessian as a NumPy array (shape
    `(3N, 3N)`, i.e. `9 × 9` for water) if you want to diagonalise or post-process
    it yourself.
  - `mol.get_results()` — the JSON-friendly summary dictionary (atoms, coords,
    `energy`, `hess`, …), matching the `<project>.json` file on disk.

The IR and Raman intensities in the log are exactly the stick spectrum: plot
intensity versus frequency to get the simulated IR/Raman spectrum for comparison
with experiment.

## Manual

- Hessian / vibrational-analysis workflow (analytical vs numerical, frequencies,
  IR/Raman, thermochemistry):
  <https://open-quantum-platform.github.io/openqp-docs/workflows/hessian/>
- `[hess]` keyword reference (`type`, `state`, `clean`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/hess/>
- `[scf]` keyword reference (`type`, `multiplicity` — choosing the reference):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- Running OpenQP from Python (the `job.workflow.hessian(...)` idiom):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
```