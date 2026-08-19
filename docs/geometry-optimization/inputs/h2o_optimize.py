"""Ground-state geometry optimization of water with the native OpenQP optimizer.

Equivalent to inputs/h2o_optimize.oqp: BHHLYP/6-31G* optimized on the native
`lib=oqp` backend with its default coordinate and trust-radius settings.
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_optimize", silent=1)

# Same starting geometry as the .oqp (Angstrom); charge 0, closed-shell singlet.
job.molecule(
    """
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
""",
    charge=0,
    multiplicity=1,
)

# Quantum theory: BHHLYP Kohn-Sham reference with the 6-31g* basis,
# matching the .oqp deck's "rks/bhhlyp/6-31g*" route.
job.theory.dft(functional="bhhlyp", basis="6-31g*", reference="rhf")

# Workflow: geometry optimization on the native optimizer, matching the
# bare `opt` of the .oqp deck.  istate=0 is the HF/DFT ground state; the
# backend (lib="oqp"), coordinate system (coordsys="auto"), trust radius
# (0.2) and cycle cap (maxit=30) are the defaults and are not written.
# Pass e.g. coordsys="tric", trust=0.1 to change one -- they are routed to
# the [oqp] backend section automatically.
job.workflow.optimize(istate=0)

mol = job.run()

# The final SCF energy is the energy at the optimized minimum.
print("Optimized SCF energy:", mol.get_scf_energy())
print("Optimized geometry (Bohr):", mol.get_system())
print(mol.get_results())
