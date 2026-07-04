"""MRSF-TDDFT excited states + gradient, BHHLYP/6-31G* on water.

Python-API equivalent of h2o_mrsf.inp: 3 MRSF roots (S0, S1, S2) and the
analytic gradient of state 2 (S1).
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_mrsf", silent=1)

# Water geometry (Angstrom), charge 0. The default multiplicity is fine here;
# job.theory.mrsf() sets the triplet ROHF *reference* multiplicity itself.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
)

# MRSF-TDDFT theory: BHHLYP/6-31G*, 3 response roots.
# Under the hood this sets [input] method=tdhf, [scf] type=rohf multiplicity=3,
# and [tdhf] type=mrsf nstate=3.
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)

# Gradient of state 2 (S1). state=... maps to [properties] grad and selects
# runtype=grad. For MRSF, state 1 is the lowest (S0-like) root, state 2 is S1.
job.workflow.gradient(state=2)

mol = job.run()

results = mol.get_results()
print("Reference (ROHF) energy:", results["energy"])
print("MRSF state total energies (S0, S1, S2):", results["td_energies"])
print("Gradient of S1 (Hartree/Bohr):", mol.get_grad())
