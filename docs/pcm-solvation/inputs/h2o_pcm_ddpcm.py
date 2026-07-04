"""Water RHF energy in a ddPCM dielectric continuum (energy-only PCM/ddX).

Python-API equivalent of h2o_pcm_ddpcm.inp. Run with:  python h2o_pcm_ddpcm.py
"""

from oqp.openqp import OpenQP

job = OpenQP("h2o_pcm", silent=1)

# Same water geometry (Angstrom) as the input deck.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# RHF reference SCF. (For the open-shell OH variant: reference="rohf", multiplicity=2.)
job.theory.hf(basis="6-31g*")

# Turn on the ddX reaction field, coupled to the reference SCF energy.
job.workflow.pcm(
    enabled=True,
    backend="ddx",
    mode="reference_scf",
    model="ddpcm",
    epsilon=78.3553,
)

mol = job.run()

# The solvated SCF energy already includes the PCM reaction-field term.
print("Solvated SCF energy (Hartree):", mol.get_scf_energy())
