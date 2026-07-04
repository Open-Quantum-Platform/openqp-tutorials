"""RHF/STO-3G GIAO NMR nuclear magnetic shielding of water.

Equivalent of inputs/h2o_giao-nmr.inp using the OpenQP Python API.
Run with:  python h2o_giao-nmr.py
"""

from oqp.openqp import OpenQP

# Same water geometry as the input deck (Angstrom).
water = """
8   0.000000000   0.000000000  -0.041061554
1  -0.533194329   0.533194329  -0.614469223
1   0.533194329  -0.533194329  -0.614469223
"""

job = OpenQP("h2o_giao_nmr", silent=1)

# Molecular identity: closed-shell singlet, neutral.
job.molecule(water, charge=0, multiplicity=1)

# Quantum theory: Hartree-Fock reference SCF in the STO-3G basis (RHF by default).
job.theory.hf(basis="sto-3g")

# Workflow: NMR shielding with the gauge-including-atomic-orbital (GIAO) gauge.
# This maps onto [properties] scf_prop=nmr, nmr_gauge=giao.
job.workflow.nmr(gauge="giao")

mol = job.run()

results = mol.get_results()
print("SCF energy (Hartree):", results["energy"])
# Per-atom isotropic shielding is the last entry of each nmr_shielding row (ppm).
for atom, row in zip(results["atoms"], results["nmr_shielding"]):
    print(f"atom Z={atom:>2}  isotropic shielding = {row[-1]:10.4f} ppm")
