"""Standalone MP2 energy for water (6-31G), matching h2o_mp2_6-31g.inp.

Run with:  python h2o_mp2_6-31g.py

Switch `variant` to "scs-mp2" or "sos-mp2" for spin-scaled MP2, and switch
`reference` to "rhf" / "uhf" / "rohf" to change the HF reference.
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_mp2", silent=1)

# Same geometry as the .inp (Angstrom). "water" is a built-in name, but we
# spell out the coordinates so the two styles produce the identical number.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# MP2 on a UHF reference; variant="mp2" is conventional (unscaled) MP2.
# conv is forwarded to the [scf] section as a tight SCF convergence.
job.theory.mp2(basis="6-31g", reference="uhf", variant="mp2", conv=1.0e-10)

mol = job.run()

results = mol.get_results()
print("HF reference energy:", mol.get_scf_energy())
print("MP2 total energy:   ", results["energy"])
