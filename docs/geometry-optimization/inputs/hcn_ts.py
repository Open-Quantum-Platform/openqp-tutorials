"""Transition-state search for HCN -> HNC isomerization via geomeTRIC.

Equivalent to inputs/hcn_ts.inp: a ground-state (method=hf) TS optimization on
the geomeTRIC backend. runtype=ts is selected by job.workflow.ts(...).
"""

from oqp.openqp import OpenQP

job = OpenQP("hcn_ts", silent=1)

# A bent starting guess near the HCN <-> HNC saddle point (Angstrom).
job.molecule(
    """
C   0.0000000000   0.0000000000   0.0000000000
N   0.0000000000   0.0000000000   1.1700000000
H  -1.1000000000   0.0000000000   0.0000000000
""",
    charge=0,
    multiplicity=1,
)

# Hartree-Fock reference with a small 3-21g basis (fast).
job.theory.hf(basis="3-21g")

# Workflow: transition-state search (runtype=ts) on the geomeTRIC backend.
#   lib="geometric"           -> [optimize] lib=geometric
#   istate=0                  -> [optimize] istate=0 (ground state)
#   coordsys/trust/tmax/...   -> routed to the [geometric] backend section
job.workflow.ts(
    lib="geometric",
    istate=0,
    maxit=50,
    coordsys="dlc",
    trust=0.05,
    tmax=0.1,
    hessian="never",
    convergence_set="GAU",
    prefix="hcn_ts",
)

mol = job.run()

# Energy at the converged saddle point.
print("TS SCF energy:", mol.get_scf_energy())
print(mol.get_results())
