"""SCF convergence with the SOSCF converger — Python-API equivalent of
h2o_scf_soscf.inp.

job.theory.dft(...) forwards any extra keyword straight into the [scf] section,
so converger_type / conv / maxit land exactly where the input-file deck puts them.
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_scf_soscf", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# DFT/PBE reference; extra kwargs are written to [scf].
# converger_type="soscf" == [scf] converger_type=soscf
job.theory.dft(
    functional="pbe",
    basis="6-31g",
    reference="rhf",       # [scf] type=rhf
    converger_type="soscf",  # second-order SCF instead of DIIS
    maxit=60,                # [scf] maxit
    conv=1.0e-6,             # [scf] conv
)

mol = job.run()

print("SCF energy:", mol.get_scf_energy())
print("Results:", mol.get_results())
