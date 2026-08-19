"""Minimum-energy conical intersection (MECI) of ethylene with MRSF-TDDFT.

Python-API equivalent of c2h4_mrsf_meci.oqp. It optimizes the twisted-pyramidal
S0/S1 crossing seam of ethylene using the native OpenQP optimizer and its
default augmented-Lagrangian crossing search.

Run with:  python c2h4_mrsf_meci.py
"""

from oqp.openqp import OpenQP

# Twisted ethylene starting geometry (Angstrom), the same coordinates as the
# input-file version. Near the S0/S1 conical intersection the CC bond is
# stretched and one CH2 is pyramidalized/twisted out of plane.
geometry = """
C  -1.6699351346837055   0.1537249235528157  -1.5459803491111643
C  -1.8079415266835852  -0.0386075716896284  -0.1602069788110266
H  -2.6609567768367581   0.2572290722092156  -2.0290359598415040
H  -1.2898503996116444  -0.7568524635289917  -2.0470428696820342
H  -1.3096398768036397   0.6557118321425524   0.5396052278505126
H  -2.3820842951209360  -0.7983813277099963   0.4308517619153288
"""

job = OpenQP("c2h4_mrsf_meci", silent=1)

# Neutral, closed-shell target states built on a triplet ROHF reference.
job.molecule(geometry, charge=0)

# MRSF-TDDFT on a BHHLYP/6-31G* triplet reference, solving 5 MRSF roots.
# State 1 is the (multiconfigurational) ground state S0, state 2 is S1.
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=5)

# Crossing-point search between state 1 (S0) and state 2 (S1).  The native
# optimizer, its default augmented-Lagrangian algorithm, and the default
# convergence thresholds are used; pass e.g. meci_search="penalty" or
# energy_gap=2e-3 only to change one of them.
job.workflow.meci(istate=1, jstate=2)

mol = job.run()

results = mol.get_results()
print("MRSF state energies at the optimized geometry:", results["td_energies"])
