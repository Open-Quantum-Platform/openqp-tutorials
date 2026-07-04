# ESPF QM/MM embedding: single-point energy and ground-state MD

Most chemistry does not happen in vacuum — a chromophore sits in solvent, a
substrate sits in an enzyme pocket, a defect sits in a crystal. **QM/MM** lets
you keep the interesting part quantum-mechanical (QM) while describing the
surroundings with a cheap classical force field (MM). This tutorial is the
entry point to QM/MM in OpenQP: it runs a **single-point embedded energy** and a
short **ground-state QM/MM molecular-dynamics** trajectory on a tiny water dimer,
where one water is QM and the other is MM. Once these two decks make sense, the
excited-state / nonadiabatic QM/MM work in the
[SOC-NAMD-QMMM tutorial](../soc-namd-qmmm/index.md) is just the same `[qmmm]`
block bolted onto a dynamics run.

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

The two decks in [`inputs/`](https://github.com/Open-Quantum-Platform/openqp-tutorials/tree/main/docs/qmmm-embedding/inputs) share one system — the water dimer in
[`water_dimer.pdb`](inputs/water_dimer.pdb), with the MM water parameterized by
the TIP3P force field in [`tip3p.xml`](inputs/tip3p.xml) — and differ only in
what they *do* with it: a single-point energy versus a short MD run.

## Input-file style

### Single-point energy — [`water_dimer_qmmm_energy.inp`](inputs/water_dimer_qmmm_energy.inp)

The QM region is the **first** water (PDB atoms 0, 1, 2); the **second** water is
MM and polarizes the QM density through ESPF. Annotated:

```ini
[input]
system     = water_dimer.pdb 0 1 2   # PDB path + 0-based indices of the QM atoms
charge     = 0
runtype    = energy                  # single-point energy (+ requested properties)
method     = hf                      # HF/DFT driver (functional below picks DFT)
functional = bhhlyp                  # half-and-half hybrid; empty => plain HF
basis      = 6-31g*
qmmm_flag  = True                    # turn on QM/MM; without it [qmmm] is ignored

[scf]
type         = rhf                   # closed-shell reference for the QM region
multiplicity = 1

[qmmm]
forcefield_files = tip3p.xml         # OpenMM force field for the MM water
cutoff           = NoCutoff          # isolated (non-periodic) cluster
embedding        = electrostatic     # full ESPF electrostatic embedding
```

Key points, section by section:

- **`system`** does double duty here. Because the first token ends in `.pdb`,
  OpenQP reads the geometry *from the PDB* and interprets the trailing numbers
  (`0 1 2`) as the **0-based indices of the QM atoms**. Everything else in the
  PDB is MM. (Ranges like `0-2` work too.)
- **`method=hf` with a non-empty `functional`** makes this a DFT (Kohn-Sham) run —
  `bhhlyp` is the half-and-half hybrid. Leave `functional` empty for plain
  Hartree-Fock.
- **`qmmm_flag=True`** is the master switch: it turns on the ESPF driver and makes
  OpenQP read `[qmmm]`. Without it the `[qmmm]` section is silently ignored and
  you get a gas-phase calculation of the QM atoms only.
- In **`[qmmm]`**, `forcefield_files` parameterizes the MM atoms (the QM atoms'
  electrostatics come from ESPF, not from fixed MM charges); `cutoff=NoCutoff`
  says isolated cluster; `embedding=electrostatic` selects the full ESPF
  electrostatic coupling (as opposed to a cheaper mechanical embedding).

### Ground-state MD — [`water_dimer_qmmm_md.inp`](inputs/water_dimer_qmmm_md.inp)

Same dimer, now **propagated in time**: the QM water moves under the embedded
QM/ESPF force, the MM water under the force field, coupled through ESPF. This is
the ground-state QM/MM MD path (the `QMMM_MD` driver, `runtype=md`). Annotated:

```ini
[input]
runtype    = md                      # ground-state QM/MM MD (OpenMM integrator)
qmmm_flag  = True
method     = hf
functional = bhhlyp
basis      = 6-31g
charge     = 0

[scf]
type         = rhf
multiplicity = 1

[qmmm]
pdb_file         = water_dimer.pdb   # full QM+MM system (MD path reads it here)
forcefield_files = tip3p.xml
qm_atoms         = 0-2               # 0-based indices of the QM water
cutoff           = NoCutoff          # isolated cluster (use PME for a periodic box)
embedding        = electrostatic     # full ESPF electrostatic embedding
n_steps          = 5                 # number of MD steps (raise for a real run)
timestep         = 0.5               # fs
ensemble         = nve               # microcanonical (Verlet); nvt/npt = Langevin
temperature      = 300.0             # K (initial velocities / thermostat target)
```

What changed versus the energy deck:

- **`runtype=md`** selects the ground-state QM/MM MD driver. (This path is *not*
  part of `openqp --run_tests all`; it needs OpenMM.)
- The QM region is now specified **inside `[qmmm]`**, not on the `system` line:
  `pdb_file` holds the full QM+MM system and `qm_atoms=0-2` carves out the QM
  water. There is no `system=` line — the MD path reads geometry from `pdb_file`.
- The extra `[qmmm]` keys are the **integrator controls**, read by the MD driver:
  `n_steps` (how many steps), `timestep` (fs), `ensemble` (`nve` microcanonical
  Verlet, or `nvt` / `npt` Langevin), and `temperature` (initial-velocity /
  thermostat target in K). Five steps at 0.5 fs is a smoke test — raise `n_steps`
  for a real trajectory.

> For **excited-state** or **nonadiabatic** QM/MM dynamics you do *not* use
> `runtype=md`; you use `runtype=namd` with an MRSF-TDDFT theory, as in the
> [SOC-NAMD-QMMM tutorial](../soc-namd-qmmm/index.md).

## Python style

The same two calculations with the compact `OpenQP` scripting interface.
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
same PDB-plus-indices `system` value as the `.inp` deck, so the two styles build
the identical calculation.

### Ground-state MD — [`water_dimer_qmmm_md.py`](inputs/water_dimer_qmmm_md.py)

```python
from oqp.openqp import OpenQP

job = OpenQP("water_dimer_qmmm_md", silent=1)

job.theory.dft(functional="bhhlyp", basis="6-31g", reference="rhf", multiplicity=1)

# On the MD path the QM region is taken from [qmmm]: pdb_file holds the full
# QM+MM system and qm_atoms selects the QM water.
job.qmmm(
    pdb_file="water_dimer.pdb",
    forcefield=["tip3p.xml"],
    qm_atoms="0-2",
    cutoff="NoCutoff",
    embedding="electrostatic",
)

mol = job.run(run_type="md")        # ground-state QM/MM MD
```

One honest caveat about the MD script: the integrator controls `n_steps`,
`timestep`, and `ensemble` are read by the `QMMM_MD` driver **directly from the
`[qmmm]` section of the input file** — they are not part of the strict Python-API
schema. So the Python script above sets up the QM/MM system and selects
`runtype=md`, but the trajectory length / step / ensemble fall back to the driver
defaults (`n_steps=1000`, `timestep=1.0` fs, `ensemble=nve`). To pin them to the
tutorial values (`n_steps=5`, `timestep=0.5`, `ensemble=nve`), run the `.inp`
with the CLI — that is the recommended entry point for MD.

## Run it

Run from the `inputs/` folder so the PDB and force-field files resolve. Either
style works for the energy deck:

```bash
cd qmmm-embedding/inputs

# single-point energy
openqp water_dimer_qmmm_energy.inp        # input-file style
python water_dimer_qmmm_energy.py         # Python-API style

# ground-state MD (CLI is the recommended entry point; see caveat above)
openqp water_dimer_qmmm_md.inp
python water_dimer_qmmm_md.py
```

Both need OpenQP (`pip install openqp`) **and** OpenMM (`pip install openmm`);
without OpenMM these decks are SKIPPED. Check both import with
`python -c "import oqp, openmm; print('ok')"`.

## Reading the output

**Single-point energy.** The number you want is the **embedded QM/MM SCF energy** —
the QM energy computed in the field of the MM charges (so it already includes the
QM–MM electrostatic coupling), not a gas-phase value.

- In the **log file** (`water_dimer_qmmm_energy.log`) look for the converged SCF
  total energy.
- From **Python**, `mol.get_scf_energy()` returns it directly. `mol.get_results()`
  gives the full results dict (matching the `<project>.json` file).

A quick sanity check that embedding is actually on: rerun with `qmmm_flag=False`
(gas-phase first water) and confirm the SCF energy shifts — the difference is the
polarization + interaction the MM water induces.

**Ground-state MD.** The `QMMM_MD` driver writes a trajectory and an energy log
as it steps:

- a **trajectory** file (PDB or DCD) with the propagated coordinates;
- an energy record (`total_energy.npz` by default) holding, per step, the
  potential / kinetic / total energy and the instantaneous temperature.

For an `nve` (microcanonical) run the diagnostic to watch is **total-energy
conservation** — `E_tot` should stay flat across the steps; a visible drift means
the step is too large or the SCF is not tight enough.

## Manual

- `[qmmm]` keyword reference (`pdb_file`, `qm_atoms`, `forcefield_files`,
  `cutoff`, `embedding`, and the MD controls):
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
