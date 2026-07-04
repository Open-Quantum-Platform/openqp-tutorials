"""Ground-state geometry optimization of water with the native OpenQP optimizer.

Equivalent to inputs/h2o_optimize.inp: RHF/BHHLYP-as-HF (method=hf) with the
6-31g* basis, optimized on the native `lib=oqp` backend in TRIC coordinates.
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_optimize", silent=1)

# Same starting geometry as the .inp (Angstrom); charge 0, closed-shell singlet.
job.molecule(
    """
O  -0.0000000000   0.0000000000  -0.0410615540
H  -0.5331943294   0.5331943294  -0.6144692230
H   0.5331943294  -0.5331943294  -0.6144692230
""",
    charge=0,
    multiplicity=1,
)

# Quantum theory: Hartree-Fock reference (method=hf) with the 6-31g* basis.
job.theory.hf(basis="6-31g*")

# Workflow: geometry optimization on the native optimizer.
#   lib="oqp"      -> [optimize] lib=oqp  (native backend)
#   istate=0       -> [optimize] istate=0 (HF/DFT ground state)
#   maxit=30       -> [optimize] maxit=30
#   coordsys/trust -> routed to the [oqp] backend section automatically
job.workflow.optimize(
    lib="oqp",
    istate=0,
    maxit=30,
    coordsys="tric",
    trust=0.2,
)

mol = job.run()

# The final SCF energy is the energy at the optimized minimum.
print("Optimized SCF energy:", mol.get_scf_energy())
print("Optimized geometry (Bohr):", mol.get_system())
print(mol.get_results())
