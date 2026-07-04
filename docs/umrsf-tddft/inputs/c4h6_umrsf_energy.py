"""UMRSF-TDDFT single-point excitation energies of butadiene (C4H6).

Equivalent to c4h6_umrsf_energy.inp: BHHLYP / 6-31G*, energy only,
a UHF triplet reference with an unrestricted-MRSF (umrsf) response.

The compact job.theory.mrsf(...) helper always sets an ROHF reference and
type=mrsf, so for the UHF/umrsf variant we set the three UMRSF-specific
sections explicitly with job.update(...) -- the same [input]/[scf]/[tdhf]
keywords the input deck uses. This raw-section style is the documented way to
express calculations the model-specific helpers do not cover.
"""

from oqp.openqp import OpenQP

# Butadiene geometry (Angstrom), same atoms/order as the .inp deck.
butadiene = """
6   -1.901080641    0.114577198    0.000000000
6   -0.574900598   -0.402237631    0.000000000
6    1.901080641   -0.114577198    0.000000000
6    0.574900598    0.402237631    0.000000000
1   -2.755248674   -0.545115910    0.000000000
1   -2.071562576    1.183086840    0.000000000
1    2.755248674    0.545115910    0.000000000
1   -0.441348756   -1.477941036    0.000000000
1    2.071562576   -1.183086840    0.000000000
1    0.441348756    1.477941036    0.000000000
"""

job = OpenQP("c4h6_umrsf", silent=1)
job.molecule(butadiene, charge=0)

# UMRSF-TDDFT = method=tdhf + UHF triplet reference + [tdhf] type=umrsf.
job.update({
    "input": {"method": "tdhf", "functional": "bhhlyp",
              "basis": "6-31g*", "runtype": "energy"},
    "scf":   {"type": "uhf", "multiplicity": 3},          # UHF high-spin reference
    "tdhf":  {"type": "umrsf", "nstate": 5, "multiplicity": 1},  # 5 singlet targets
})

# Serial run for an ordinary (non-MPI) Python invocation.
job.control(usempi=False)

mol = job.run()
results = mol.get_results()

print("Reference (SCF) energy:", results["energy"])
print("UMRSF state energies (Hartree, relative to reference):")
for i, e in enumerate(results["td_energies"]):
    print(f"  state {i}: {e:+.6f}")
