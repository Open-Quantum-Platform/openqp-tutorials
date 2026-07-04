"""HF and DFT energies and gradients on water via the OpenQP Python API.

The first block reproduces inputs/h2o_rhf_energy.inp exactly (RHF/6-31G*
single-point energy). The remaining blocks show the ROHF/UHF references, the
DFT (Kohn-Sham) helper, and how to ask for a gradient instead of an energy.
Run with:  python h2o_rhf_energy.py
"""

from oqp.openqp import OpenQP

# Same water geometry (Angstrom) used by the .inp files.
WATER = """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""

# --- RHF/6-31G* energy: the exact equivalent of h2o_rhf_energy.inp ----------
job = OpenQP("h2o_rhf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.hf(reference="rhf", basis="6-31g*")   # closed-shell HF
mol = job.run()
print("RHF energy:", mol.get_scf_energy())        # -> -76.0107 Hartree

# --- ROHF triplet reference (open-shell, restricted) ------------------------
job = OpenQP("h2o_rohf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=3)
job.theory.hf(reference="rohf", basis="6-31g*")
print("ROHF energy:", job.run().get_scf_energy())

# --- UHF triplet reference (open-shell, unrestricted) -----------------------
job = OpenQP("h2o_uhf_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=3)
job.theory.hf(reference="uhf", basis="6-31g*", stability=True)
print("UHF energy:", job.run().get_scf_energy())

# --- DFT (BHHLYP) energy: job.theory.dft instead of job.theory.hf -----------
job = OpenQP("h2o_dft_energy", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.dft(functional="bhhlyp", reference="rhf", basis="6-31g*")
print("DFT (BHHLYP) energy:", job.run().get_scf_energy())

# --- RHF gradient: swap the energy workflow for job.workflow.gradient -------
# state=0 is the SCF (ground-state) gradient; it maps to [properties] grad=0.
job = OpenQP("h2o_rhf_grad", silent=1)
job.molecule(WATER, charge=0, multiplicity=1)
job.theory.hf(reference="rhf", basis="6-31g*")
job.workflow.gradient(state=0)
mol = job.run()
print("RHF energy:", mol.get_scf_energy())
print("RHF gradient (Hartree/Bohr):", mol.get_grad())
