#!/usr/bin/env python
"""ESPF QM/MM single-point energy of a water dimer, via the OpenQP Python API.

The exact same calculation as water_dimer_qmmm_energy.inp. The QM region is the
first water (PDB atoms 0,1,2); the second water is MM (TIP3P point charges) and
polarizes the QM density through the ESPF operator. Requires the optional OpenMM
backend (pip install openmm). Run from this folder so the PDB / force-field
files resolve:  python water_dimer_qmmm_energy.py
"""
from oqp.openqp import OpenQP

job = OpenQP("water_dimer_qmmm_energy", silent=1)

# QM geometry + atom selection come from the PDB: "<file>.pdb <0-based indices>".
# Here atoms 0,1,2 (the first water) are QM; the rest of the PDB is MM.
job.molecule("water_dimer.pdb 0 1 2", charge=0)

# HF/DFT reference for the QM region. functional="bhhlyp" makes it DFT (the
# half-and-half hybrid); pass functional="" or drop it for plain HF.
# reference="rhf" is the closed-shell reference (multiplicity 1).
job.theory.dft(functional="bhhlyp", basis="6-31g*", reference="rhf", multiplicity=1)

# ESPF QM/MM embedding: forcefield_files parameterizes the MM water,
# cutoff=NoCutoff is an isolated (non-periodic) cluster, embedding=electrostatic
# is the full ESPF electrostatic coupling. job.qmmm sets [input] qmmm_flag=True.
job.qmmm(
    forcefield=["tip3p.xml"],
    cutoff="NoCutoff",
    embedding="electrostatic",
)

# Default runtype is energy (single-point). Run and read the SCF energy back.
mol = job.run()
print("Embedded QM/MM SCF energy:", mol.get_scf_energy())
