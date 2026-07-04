"""RHF/6-31G* single-point energy of water, with SCF population/moment analysis.

Python-API equivalent of h2o_properties.inp: it requests the electric dipole and
the Mulliken, Lowdin, and RESP atomic charges, then prints them.

Run with:  python h2o_properties.py
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_properties", silent=1)

# Same geometry (Angstrom) as the input deck.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# RHF reference at 6-31G*. theory.hf() sets [input] method=hf and [scf] type=rhf.
job.theory.hf(basis="6-31g*", reference="rhf")

# Opt-in SCF properties: [properties] scf_prop=el_mom,mulliken,lowdin,resp.
# job.settings.<section>(...) writes raw OpenQP section keywords.
job.settings.properties(scf_prop="el_mom,mulliken,lowdin,resp")

mol = job.run()
results = mol.get_results()

print("SCF energy (Ha):", results["energy"])
print("Dipole (a.u.):", results["dipole"])
print("Mulliken charges:", results["mulliken_charges"])
print("Lowdin charges:", results["lowdin_charges"])
print("RESP charges:", results["resp_charges"])
