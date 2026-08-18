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

### Single-point energy — [`water_dimer_qmmm_energy.oqp`](inputs/water_dimer_qmmm_energy.oqp)

The QM region is the **first** water (PDB atoms 0, 1, 2); the **second** water is
MM and polarizes the QM density through ESPF. Annotated:

```text
rks/bhhlyp/6-31g*                   # QM level of theory (Kohn-Sham, BHHLYP)
geom="water_dimer.pdb 0 1 2"        # PDB path + 0-based indices of the QM atoms
energy()                            # single-point energy (+ requested properties)
qmmm(forcefield_files=tip3p.xml,cutoff=NoCutoff,embedding=electrostatic)
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
  electrostatics come from ESPF, not from fixed MM charges); `cutoff=NoCutoff`
  says isolated cluster; `embedding=electrostatic` selects the full ESPF
  electrostatic coupling (as opposed to a cheaper mechanical embedding).

### Ground-state MD — [`water_dimer_qmmm_md.oqp`](inputs/water_dimer_qmmm_md.oqp)

Same dimer, now **propagated in time**: the QM water moves under the embedded
QM/ESPF force, the MM water under the force field, coupled through ESPF. This is
the ground-state QM/MM MD path (the `QMMM_MD` driver). Annotated:

```text
rks/bhhlyp/6-31g
geom="water_dimer.pdb 0-2"          # full system + the QM water's 0-based indices
md()                                # ground-state QM/MM MD (OpenMM integrator)
qmmm(forcefield_files=tip3p.xml,cutoff=NoCutoff,embedding=electrostatic,
     n_steps=5,timestep=0.5,ensemble=nve,temperature=300.0)
```

What changed versus the energy deck:

- **`md()` replaces `energy()`**, selecting the ground-state QM/MM MD driver.
  (This path is *not* part of `openqp --run-tests all`; it needs OpenMM.)
- **The QM region is still named the same way** — `geom="water_dimer.pdb 0-2"` —
  even though the legacy MD deck specified it with a separate `[qmmm] pdb_file`
  plus `qm_atoms`. One geometry line serves both drivers, and the lowering fills
  in `pdb_file` and `qm_atoms` for the MD path.
- The extra `qmmm(...)` keys are the **integrator controls**, read by the MD
  driver: `n_steps` (how many steps), `timestep` (fs), `ensemble` (`nve`
  microcanonical Verlet, or `nvt` / `npt` Langevin), and `temperature`
  (initial-velocity / thermostat target in K). Five steps at 0.5 fs is a smoke
  test — raise `n_steps` for a real trajectory.

> For **excited-state** or **nonadiabatic** QM/MM dynamics you do *not* use
> `md()`; you use `namd(...)` with an MRSF-TDDFT route, as in the
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
same PDB-plus-indices `system` value as the `.oqp` deck, so the two styles build
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
QM/MM section of the input file** — they are not part of the strict Python-API
schema. So the Python script above sets up the QM/MM system and selects the MD
run type, but the trajectory length / step / ensemble fall back to the driver
defaults (`n_steps=1000`, `timestep=1.0` fs, `ensemble=nve`). To pin them to the
tutorial values (`n_steps=5`, `timestep=0.5`, `ensemble=nve`), run the `.oqp` deck
with the CLI — that is the recommended entry point for MD.

## Run it

Run from the `inputs/` folder so the PDB and force-field files resolve. Either
style works for the energy deck:

```bash
cd qmmm-embedding/inputs

# single-point energy
openqp water_dimer_qmmm_energy.oqp        # input-file style
python water_dimer_qmmm_energy.py         # Python-API style

# ground-state MD (CLI is the recommended entry point; see caveat above)
openqp water_dimer_qmmm_md.oqp
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
