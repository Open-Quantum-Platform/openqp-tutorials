"""MRSF-TDDFT spin-orbit coupling of water via the OpenQP Python API.

Equivalent to inputs/h2o_soc.inp: a triplet-ROHF MRSF-TDDFT calculation that
forms the singlet/triplet manifold and the SOC matrix elements between them.
Run with:  python h2o_soc.py
"""

from oqp.openqp import OpenQP

# Cartesian AOs (ispher=false): 2e SOC is not available with spherical AOs.
job = OpenQP("h2o_soc", silent=1, input={"ispher": "false"})

# Same water geometry as the input deck (Angstrom).
job.molecule(
    "8  0.000000  0.000000  0.000000; "
    "1  0.772598  0.555678  0.000000; "
    "1 -0.773128  0.555678  0.000000",
    charge=0,
)

# MRSF-TDDFT theory: triplet ROHF reference + [tdhf] type=mrsf, 12 roots.
job.theory.mrsf(functional="bhhlyp", basis="6-31G(2df,p)", nstate=12)

# SOC workflow: sets runtype=soc, soc_2e, and [scf] scal_rel.
#   soc_2e=1  -> 1e Breit-Pauli + mean-field 2e SOC
#   scal_rel=2 -> DKH scalar-relativistic core (helper default is 2)
job.workflow.soc(soc_2e=1, scal_rel=2)

mol = job.run()

# SOC eigenvalues (spin-adiabatic states) in cm^-1.
print("SOC eigenvalues (cm^-1):", mol.get_soc())

# Full JSON-friendly summary: td_singlet_energies, td_triplet_energies, soc, ...
print(mol.get_results())
