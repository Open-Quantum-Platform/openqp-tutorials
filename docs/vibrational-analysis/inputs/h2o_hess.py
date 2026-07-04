"""Analytical HF/DFT Hessian of water: frequencies, normal modes, IR/Raman.

Equivalent to inputs/h2o_hess.inp. Run with:  python h2o_hess.py
"""

from oqp.openqp import OpenQP

# 1. Create the job (silent=1 keeps the console quiet; the .log still records everything).
job = OpenQP("h2o_hess", silent=1)

# 2. Molecular identity. The built-in "water" geometry matches the deck; you can
#    also pass an inline XYZ block or an .xyz path here.
job.molecule(geometry="water", charge=0, multiplicity=1)

# 3. Quantum theory: BHHLYP/6-31G* Kohn-Sham. job.theory.dft(...) sets
#    method=hf + functional=bhhlyp + [scf] type=rhf, exactly like the deck.
job.theory.dft(functional="bhhlyp", basis="6-31g*")

# 4. Workflow: the analytical ground-state Hessian. This sets runtype=hess and
#    the [hess] section. state=0 is the HF/DFT ground state; clean=True removes
#    the scratch Hessian files afterwards.
job.workflow.hessian(type="analytical", state=0, clean=True)

# 5. Run and collect results.
mol = job.run()

print("SCF energy (Hartree):", mol.get_scf_energy())

# The mass-weighted Hessian is diagonalised into harmonic frequencies, normal
# modes, thermochemistry, and IR/Raman intensities in the .log. The raw
# Cartesian Hessian matrix is available here for your own post-processing.
hessian = mol.get_hess()
print("Hessian shape:", None if hessian is None else hessian.shape)

# get_results() is the JSON-friendly summary (atoms, coords, energy, hess, ...).
print(mol.get_results().keys())
