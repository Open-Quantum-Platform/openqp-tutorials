# ESPF QM/MM embedding: a single-point embedded energy

Most chemistry does not happen in vacuum — a chromophore sits in solvent, a
substrate sits in an enzyme pocket, a defect sits in a crystal. **QM/MM** lets
you keep the interesting part quantum-mechanical (QM) while describing the
surroundings with a cheap classical force field (MM). This tutorial is the
entry point to QM/MM in OpenQP: it runs a **single-point embedded energy** on a
tiny water dimer, where one water is QM and the other is MM. Once this deck makes
sense, the excited-state / nonadiabatic QM/MM work in the
[SOC-NAMD-QMMM tutorial](../soc-namd-qmmm/index.md) is just the same `qmmm(...)`
call bolted onto a dynamics run.

## A little theory

OpenQP couples the QM and MM regions with **ESPF electrostatic embedding**
(Electrostatic-Potential-Fitted). The MM atoms carry fixed point charges from
the force field; those charges create an electrostatic potential at the QM
region, and ESPF folds that potential into the QM Hamiltonian through a small set
of fitted one-electron operators. Two things follow:

- the **MM charges polarize the QM density** — the QM wavefunction relaxes in the
  field of its environment, exactly as it would if the environment were quantum;
- the coupling has an **analytic, energy-conserving gradient**, so the same
  machinery drives geometry optimization and molecular dynamics.

`cutoff=NoCutoff` treats the system as an isolated (non-periodic) cluster; for a
solvated periodic box you switch to `cutoff=PME`, which evaluates the embedding
with particle-mesh Ewald. The MM engine underneath is [OpenMM](https://openmm.org),
so QM/MM decks need it installed (`pip install openmm`); without it they are
reported **SKIPPED** rather than run. For the full derivation and input contract
see the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

The deck in [`inputs/`](https://github.com/Open-Quantum-Platform/openqp-tutorials/tree/main/docs/qmmm-embedding/inputs) uses the water dimer in
[`water_dimer.pdb`](inputs/water_dimer.pdb), with the MM water parameterized by
the TIP3P force field in [`tip3p.xml`](inputs/tip3p.xml).

> **Ground-state QM/MM MD** used to have a second deck here. OpenQP's input
> checker currently reports `runtype=md` as *recognized but not implemented*, so
> that deck could not be run as written and has been removed rather than left as
> a trap. For QM/MM dynamics today, use the nonadiabatic `namd(...)` path in the
> [SOC-NAMD-QMMM tutorial](../soc-namd-qmmm/index.md).

## Input-file style

### Single-point energy — [`water_dimer_qmmm_energy.oqp`](inputs/water_dimer_qmmm_energy.oqp)

The QM region is the **first** water (PDB atoms 0, 1, 2); the **second** water is
MM and polarizes the QM density through ESPF. Annotated:

```text
rks/bhhlyp/6-31g* qmmm(forcefield_files="tip3p.xml")   # QM level of theory + the MM force field
geom="water_dimer.pdb 0 1 2"        # PDB path + 0-based indices of the QM atoms
```

Key points:

- **`geom` does double duty here.** Because the value ends in `.pdb` followed by
  indices, OpenQP reads the geometry *from the PDB* and interprets the trailing
  numbers (`0 1 2`) as the **0-based indices of the QM atoms**. Everything else in
  the PDB is MM. (Ranges like `0-2` work too.)
- **`rks/bhhlyp/6-31g*`** makes the QM region a Kohn-Sham DFT calculation with the
  half-and-half hybrid. Drop the functional component (`rhf/6-31g*`) for plain
  Hartree-Fock.
- **`qmmm(...)` is the master switch.** Supplying it turns on the ESPF driver *and*
  carries the MM settings; there is no separate flag to forget, and the PDB path
  from `geom` is propagated into the section for you. (A top-level `qmmm=true`
  exists for turning QM/MM on with defaults, but combining it with `qmmm(...)` is
  rejected as saying the same thing twice.)
- Inside the call, `forcefield_files` parameterizes the MM atoms (the QM atoms'
  electrostatics come from ESPF, not from fixed MM charges). The defaults —
  `cutoff=NoCutoff` (isolated cluster) and `embedding=electrostatic` (full ESPF
  electrostatic coupling, as opposed to a cheaper mechanical embedding) — are
  what this cluster needs, so they are not written.
- **No driver keyword means a single-point energy.**

## Python style

The same calculation with the compact `OpenQP` scripting interface.
`job.qmmm(...)` sets `[input] qmmm_flag=True` and fills `[qmmm]`; `job.theory.dft(...)`
sets `method`/`functional` and the `[scf]` reference.

### Single-point energy — [`water_dimer_qmmm_energy.py`](inputs/water_dimer_qmmm_energy.py)

```python
from oqp.openqp import OpenQP

job = OpenQP("water_dimer_qmmm_energy", silent=1)

# QM geometry + atom selection come from the PDB: "<file>.pdb <0-based indices>".
# Atoms 0,1,2 (the first water) are QM; the rest of the PDB is MM.
job.molecule("water_dimer.pdb 0 1 2", charge=0)

# bhhlyp -> DFT (half-and-half hybrid); reference="rhf" is closed-shell.
job.theory.dft(functional="bhhlyp", basis="6-31g*", reference="rhf", multiplicity=1)

# ESPF QM/MM embedding. job.qmmm sets [input] qmmm_flag=True.
job.qmmm(
    forcefield=["tip3p.xml"],       # alias for [qmmm] forcefield_files
    cutoff="NoCutoff",              # isolated cluster; "PME" for a periodic box
    embedding="electrostatic",
)

mol = job.run()                     # default runtype is energy
print("Embedded QM/MM SCF energy:", mol.get_scf_energy())
```

The `"water_dimer.pdb 0 1 2"` string passed to `job.molecule(...)` is the exact
same PDB-plus-indices `system` value as the `.oqp` deck, so the two styles build
the identical calculation.

## Run it

Run from the `inputs/` folder so the PDB and force-field files resolve. Either
style works:

```bash
cd qmmm-embedding/inputs
openqp water_dimer_qmmm_energy.oqp        # input-file style
python water_dimer_qmmm_energy.py         # Python-API style
```

Both need OpenQP (`pip install openqp`) **and** OpenMM (`pip install openmm`);
without OpenMM the deck is SKIPPED. Check both import with
`python -c "import oqp, openmm; print('ok')"`.

## Reading the output

The number you want is the **embedded QM/MM SCF energy** — the QM energy computed
in the field of the MM charges (so it already includes the QM–MM electrostatic
coupling), not a gas-phase value.

- In the **log file** (`water_dimer_qmmm_energy.log`) look for the converged SCF
  total energy.
- From **Python**, `mol.get_scf_energy()` returns it directly. `mol.get_results()`
  gives the full results dict (matching the `<project>.json` file).

A quick sanity check that embedding is actually on: rerun with `qmmm_flag=False`
(gas-phase first water) and confirm the SCF energy shifts — the difference is the
polarization + interaction the MM water induces.

## Manual

- `[qmmm]` keyword reference (`pdb_file`, `qm_atoms`, `forcefield_files`,
  `cutoff`, `embedding`):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/qmmm/>
- QM/MM overview and ESPF embedding (theory + input contract, PME/periodic setup,
  covalent boundaries):
  <https://open-quantum-platform.github.io/openqp-docs/>
- Running OpenQP from Python (the `job.qmmm(...)` / `job.theory.dft(...)` idiom):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>

## References

- ESPF QM/MM embedding — Huix-Rotllant & Ferré, [10.1021/acs.jctc.0c01075](https://doi.org/10.1021/acs.jctc.0c01075)
- OpenMM — Eastman et al., [10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659)
- OpenQP — [10.1021/acs.jctc.4c01117](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)
