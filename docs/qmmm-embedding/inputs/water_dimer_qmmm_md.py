#!/usr/bin/env python
"""Ground-state ESPF QM/MM molecular dynamics of a water dimer (runtype=md).

Python-API counterpart of water_dimer_qmmm_md.inp. The QM water (atoms 0-2) is
propagated under the embedded QM/ESPF force and the MM water under the force
field, coupled through ESPF electrostatics, using an OpenMM integrator. This is
the ground-state QM/MM MD path (the QMMM_MD driver); for excited-state /
nonadiabatic QM/MM dynamics use job.workflow.namd(...) with runtype=namd instead
(see the SOC-NAMD-QMMM tutorial). Requires OpenMM. Run from this folder:
    python water_dimer_qmmm_md.py

NOTE on the MD controls. n_steps / timestep / ensemble are read by the QMMM_MD
driver directly from the [qmmm] section of the *input file*; they are not part
of the strict Python-API schema. So this script sets up the QM/MM system with
job.qmmm(...) and selects runtype=md, and the trajectory length / step / ensemble
fall back to the driver defaults (n_steps=1000, timestep=1.0 fs, ensemble=nve).
To pin them to the values in water_dimer_qmmm_md.inp (n_steps=5, timestep=0.5,
ensemble=nve), run that .inp with the CLI:  openqp water_dimer_qmmm_md.inp
"""
from oqp.openqp import OpenQP

job = OpenQP("water_dimer_qmmm_md", silent=1)

# HF/DFT reference for the QM region (bhhlyp -> DFT), closed-shell RHF.
job.theory.dft(functional="bhhlyp", basis="6-31g", reference="rhf", multiplicity=1)

# ESPF QM/MM. On the MD path the QM region and geometry are taken from the
# [qmmm] section: pdb_file holds the full QM+MM system and qm_atoms selects the
# QM water. job.qmmm sets [input] qmmm_flag=True and fills [qmmm].
job.qmmm(
    pdb_file="water_dimer.pdb",     # full QM+MM system
    forcefield=["tip3p.xml"],
    qm_atoms="0-2",                 # 0-based indices of the QM water
    cutoff="NoCutoff",              # isolated cluster; use PME for a periodic box
    embedding="electrostatic",      # full ESPF electrostatic embedding
)

# Select the ground-state QM/MM MD workflow (runtype=md) and propagate.
mol = job.run(run_type="md")
