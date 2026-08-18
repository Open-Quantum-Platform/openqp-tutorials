#!/usr/bin/env python
"""SOC-NAMD-QMMM of formaldehyde in water, via the OpenQP Python API.

The exact same calculation as h2co-water_soc-namd-qmmm.oqp, built with the
compact `OpenQP` scripting interface. Run from this folder (so the PDB /
force-field files resolve):  python h2co-water_soc-namd-qmmm.py
"""
from oqp.openqp import OpenQP

job = OpenQP("h2co-water_soc-namd-qmmm", silent=1)

# QM subsystem geometry (the 4 formaldehyde atoms selected by qm_atoms below).
job.molecule(
    """
    C   0.000000   0.000000   0.000000
    O   0.000000   0.000000   1.203000
    H   0.000000   0.943000  -0.589000
    H   0.000000  -0.943000  -0.589000
    """,
    charge=0,
)

# MRSF-TDDFT on a triplet (multiplicity=3) ROHF reference; 2 singlet roots.
# job.theory.mrsf sets [scf] type=rohf multiplicity=3 and [tdhf] type=mrsf.
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=2, multiplicity=3)

# ESPF QM/MM: the PDB holds all atoms, qm_atoms picks the QM region, the rest is
# MM. `forcefield` is the alias for [qmmm] forcefield_files.
job.qmmm(
    pdb_file="formaldehyde_water.pdb",
    forcefield=["formaldehyde.xml", "tip3p.xml"],
    qm_atoms="0-3",
    cutoff="NoCutoff",          # isolated cluster; use "PME" for a periodic box
    embedding="electrostatic",
)

# Surface hopping with spin-orbit coupling ON (soc=True -> intersystem crossing).
# job.workflow.namd selects runtype=namd and fills the [md] section.
job.workflow.namd(
    soc=True,
    nstep=1,                    # raise (e.g. 500-2000) for a real trajectory
    dt=0.25,
    active=5,                   # initially populated spin-adiabatic state
    substep=50,
    init_temp=300,
    velocity="maxwell",
    seed=3,
    decoherence="edc",
    trivial=True,
    thrshe=0.1,
)

mol = job.run()
