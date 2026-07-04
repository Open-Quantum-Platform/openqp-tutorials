"""TDDFT / 6-31G* / B3LYP5 vertical excitation energies of H2O.

Equivalent to inputs/h2o_tddft_energy.inp, driven through the OpenQP Python API.
job.theory.tddft(...) sets method=tdhf with a functional (TDDFT) and an RHF
reference; nstate is the number of excited-state roots to solve.

For plain TDHF, swap the theory line for:
    job.theory.tdhf(reference="rhf", basis="6-31g*", nstate=3)
(no functional -> Hartree-Fock linear response).
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_tddft", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# TDDFT: RHF reference + B3LYP5 functional, 3 excited-state roots.
job.theory.tddft(functional="b3lyp5", basis="6-31g*", nstate=3)

mol = job.run()
results = mol.get_results()

print("Ground-state (RHF) energy [Ha]:", results["energy"])
print("Excitation energies [Ha]     :", results["td_energies"])
