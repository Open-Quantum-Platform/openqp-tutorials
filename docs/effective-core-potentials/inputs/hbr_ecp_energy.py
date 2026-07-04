"""BHHLYP/aug-cc-pVDZ-PP energy of HBr with an ECP on Br (Python API).

Equivalent to hbr_ecp_energy.inp. The "-PP" basis carries an effective core
potential for bromine's core electrons; hydrogen keeps an all-electron basis.
The per-atom basis list is the documented ECP idiom: OpenQP reads the ECP that
ships with the "-PP" set and builds it through libecpint.
"""

from oqp.openqp import OpenQP

job = OpenQP("hbr_ecp_energy", silent=1)

# Geometry in Angstrom, Br then H (bond length 1.407611 Angstrom).
job.molecule(
    """
Br  0.000000000  0.000000000  0.000000000
H   0.000000000  0.000000000  1.407611
""",
    charge=0,
    multiplicity=1,
)

# Closed-shell BHHLYP (DFT uses the .dft helper; method=hf under the hood).
job.theory.dft(functional="bhhlyp")

# Per-atom basis, in atom order: ECP-carrying set on Br, all-electron on H.
job.settings.basis(["aug-cc-pVDZ-PP", "aug-cc-pVDZ"])

mol = job.run()

print("SCF energy:", mol.get_scf_energy())
print("Results:", mol.get_results())
