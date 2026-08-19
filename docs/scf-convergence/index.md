# SCF convergence, initial guesses, and TRAH: getting the reference to converge

Every wavefunction and DFT calculation in OpenQP begins by solving the
**self-consistent field (SCF)** equations: you guess a set of molecular orbitals,
build the Fock matrix they imply, diagonalize it to get *new* orbitals, and repeat
until the orbitals stop changing. Most closed-shell molecules converge in a dozen
iterations without any thought. But stretched bonds, transition-metal complexes,
diradicals, and near-degenerate states can make the default extrapolation stall,
oscillate, or slide into the wrong solution. This tutorial shows the two knobs that
fix almost all of those cases — the **initial guess** (`guess(type=...)`) and the
**converger** (`scf(converger_type=...)`, including **SOSCF** and **TRAH**) — on a
plain water/PBE calculation you can run in a second.

## A little theory

The SCF cycle is a fixed-point iteration, and a naive "plug the new orbitals
straight back in" scheme frequently diverges. OpenQP gives you a ladder of
convergers, each doing more work per iteration but converging more robustly:

- **DIIS** (the default) extrapolates the next Fock matrix from a history of
  previous ones. Cheap and excellent for well-behaved, closed-shell cases.
- **SOSCF** (second-order SCF) uses curvature — an approximate orbital-rotation
  Hessian — to take Newton-like steps once you are near the solution, so it
  finishes the tail of a hard convergence faster than DIIS alone. It is the
  converger used in this tutorial.
- **TRAH** (trust-region augmented-Hessian SCF) is the heavy artillery: a genuine
  second-order method wrapped in a trust region so the Newton step is only taken as
  far as it can be trusted. It converges cases where DIIS and SOSCF both fail —
  strongly correlated references, tricky open shells — at a higher per-iteration
  cost.

The other half of the battle is where you *start*. A better **initial guess** puts
you inside the basin of the right solution, which cuts iterations and avoids
converging to an excited or symmetry-broken state. OpenQP's default is the
extended-**Hückel** guess (`guess(type=huckel)`); it also offers `hcore`, `sap`, `minao`,
`modhuckel`, and JSON-restart guesses. For the full contract behind each option,
see the [`[scf]`](https://open-quantum-platform.github.io/openqp-docs/keywords/scf/)
and [`[guess]`](https://open-quantum-platform.github.io/openqp-docs/keywords/guess/)
keyword pages.

## Input-file style

The runnable deck is [`inputs/h2o_scf_soscf.oqp`](inputs/h2o_scf_soscf.oqp) — water
in the 6-31G basis, a spin-restricted PBE reference, converged with the **SOSCF**
converger starting from a **Hückel** guess. Annotated:

```text
rks/pbe/6-31g                              # restricted Kohn-Sham, PBE, 6-31G
guess(save_mol=true)                       # keep the converged orbitals for restarts
scf(converger_type=soscf,maxit=60)         # how the SCF gets there
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Key points:

- **`rks/pbe/6-31g` is Kohn-Sham DFT, not bare HF.** The middle component is the
  functional; `rhf/6-31g` — no functional at all — would give a Hartree-Fock
  reference. `uks` and `roks` are the unrestricted and restricted-open-shell
  Kohn-Sham models.
- **`guess(...)` chooses where the SCF starts.** The extended-Huckel guess is the
  default and is not written; `guess(type=hcore)`, `sap`, or `minao` changes only
  the starting orbitals, not the final answer, but can change how many iterations
  (and whether) you converge. `save_mol=true` writes the converged orbitals so a
  later run can restart from them.
- **`scf(converger_type=...)` chooses how it gets there.** This deck uses `soscf`.
  Change one option to switch strategy:

  ```text
  scf(maxit=60)                            # default first-order DIIS extrapolation
  ```

  ```text
  scf(converger_type=trah,maxit=60)        # trust-region augmented-Hessian: for hard, stalling cases
  ```

  All three converge this easy water case to the same energy; the difference shows
  up on difficult references, where `diis` may stall and `trah` still converges.
- **`maxit` and `conv` are the stopping rules.** Raise `maxit` (here 60, above the
  default of 30) for slow cases; add `conv=1e-8` when you need high-precision
  energies or gradients (the default is `1e-6`).
- **`dftgrid(...)` sets the XC integration grid.** The default (96 radial × 302
  angular points, SG2 pruning) is fine; a coarse grid can itself cause the SCF
  energy to jitter and fail to converge, so grid quality and SCF convergence are
  linked.

Note the shape of the deck: the route line says *what physics*, and everything
that follows is a named knob on one legacy section. Convergence controls are
exactly the kind of implementation detail the concise surface leaves in exact
section calls rather than promoting to the route.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_scf_soscf.py`](inputs/h2o_scf_soscf.py). `job.theory.dft(...)` sets the
`[input] functional`, the basis, and the `[scf]` reference in one call; any **extra
keyword** (here `converger_type`, `maxit`, `conv`) is forwarded straight into the
`[scf]` section, so it lands exactly where the input-file deck puts it.

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_scf_soscf", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# DFT/PBE reference; extra kwargs are written to [scf].
# converger_type="soscf" == [scf] converger_type=soscf
job.theory.dft(
    functional="pbe",
    basis="6-31g",
    reference="rhf",         # [scf] type=rhf
    converger_type="soscf",  # second-order SCF instead of DIIS
    maxit=60,                # [scf] maxit
    conv=1.0e-6,             # [scf] conv
)

mol = job.run()

print("SCF energy:", mol.get_scf_energy())
print("Results:", mol.get_results())
```

To try a different converger, change only the forwarded keyword:

```python
job.theory.dft(functional="pbe", basis="6-31g", reference="rhf", converger_type="diis")   # DIIS
job.theory.dft(functional="pbe", basis="6-31g", reference="rhf", converger_type="trah")   # TRAH
```

The full script is [`inputs/h2o_scf_soscf.py`](inputs/h2o_scf_soscf.py).

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd scf-convergence/inputs
openqp h2o_scf_soscf.oqp
```

Python style:

```bash
cd scf-convergence/inputs
python h2o_scf_soscf.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same SCF energy.

## Reading the output

This is a single-point SCF run, so the number you want is the **converged SCF
(Kohn-Sham) total energy**.

- In the **log file** (`<project>.log`) watch the **iteration table**: each row
  prints the energy and the density/gradient change, and you want to see that
  change shrink monotonically past `conv=1.0e-6` and the run report **converged**
  within `maxit` iterations. If it hits `maxit` without converging, that is the
  signal to switch `converger_type` (try `trah`) or improve the `guess(...)`.
- From **Python**, `mol.get_scf_energy()` returns the converged SCF total energy,
  and `mol.get_results()` returns the results dictionary (the `energy` field
  matches what is written to `<project>.json`).
- Because `[guess] save_mol=true`, the converged orbitals are also saved, so a
  follow-up calculation can restart from them (via a `json`/`auto` guess) instead of
  re-converging from Hückel.

The physical energy is identical no matter which converger you pick — DIIS, SOSCF,
and TRAH only differ in *how reliably and how fast* they reach it. The point of this
tutorial is to know which knob to turn when the default DIIS run refuses to converge.

## Manual

- `[scf]` keyword reference (`converger_type` = `diis` | `soscf` | `trah`, plus
  `type`, `multiplicity`, `maxit`, `conv`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- `[guess]` keyword reference (`type` = `huckel` | `hcore` | `sap` | `minao` | … ,
  and `save_mol` restart): <https://open-quantum-platform.github.io/openqp-docs/keywords/guess/>
- `[dftgrid]` keyword reference (`rad_npts`, `ang_npts`, `pruned`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/dftgrid/>
- HF / DFT workflow overview:
  <https://open-quantum-platform.github.io/openqp-docs/workflows/hf-dft/>
