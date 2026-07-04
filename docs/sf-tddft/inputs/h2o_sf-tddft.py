"""Spin-flip TDDFT on H2O with the OpenQP Python API.

Equivalent to inputs/h2o_sf-tddft.inp: a high-spin (triplet) ROHF reference is
spin-flipped down to the singlet manifold with BHHLYP/6-31G*, then the gradient
of the lowest SF root is computed.
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_sf_tddft", silent=1)

# Water geometry (Angstrom), same as the input deck.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
)

# Spin-flip TDDFT: the helper sets [tdhf] type=sf and defaults to the
# high-spin reference (reference="rohf", multiplicity=3) that SF requires.
job.theory.sf_tddft(functional="bhhlyp", basis="6-31g*", nstate=3)

# Gradient of SF root 3. For SF/MRSF, state 1 is the lowest root (the
# multiconfigurational ground state); this maps to [properties] grad=3.
job.workflow.gradient(state=3)

mol = job.run()

results = mol.get_results()
print("Reference (ROHF) energy:", results["energy"])
print("SF-TDDFT root energies:", results["td_energies"])
print("Gradient (root 3):", mol.get_grad())
