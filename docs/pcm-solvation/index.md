# PCM/ddX: implicit solvation for RHF/ROHF energies

Most chemistry happens in solution, but explicitly surrounding your molecule
with hundreds of solvent molecules is expensive and needs averaging over many
configurations. **Implicit (continuum) solvation** sidesteps that: it replaces
the solvent with a structureless **polarizable dielectric** that surrounds the
solute cavity and reacts back on the solute's charge density. This tutorial runs
OpenQP's production continuum path — an **energy-only reference-SCF PCM/ddX**
calculation — on a closed-shell (RHF) water molecule, and shows the equivalent
open-shell (ROHF) radical. Reach for it whenever you want a solvated SCF energy
(solvation shifts, relative energies in solution, a cheap solvent correction on
top of a gas-phase geometry) without the cost of explicit solvent.

## A little theory

The **Polarizable Continuum Model (PCM)** carves a cavity around the solute,
fills the space outside with a dielectric of permittivity ε (78.36 for water),
and solves for the **apparent surface charge** the solute's field induces on the
cavity boundary. That induced charge produces a **reaction-field potential**
which is added back into the SCF, so the solute and the continuum polarize each
other self-consistently. OpenQP builds this through **ddX**, a
*domain-decomposition* solver: the cavity is split into per-atom spherical
domains and the boundary problem is solved locally and stitched together, which
scales far better than a global surface mesh. ddX offers three flavors —
**ddCOSMO** (conductor limit), **ddPCM** (finite-ε dielectric, used here), and
**ddLPB** (linearized Poisson-Boltzmann, for ionic solvents). In the current
production path the reaction field couples to the **reference SCF energy** only;
excited-state PCM, PCM gradients, and PCM geometry optimization are future
extensions. For the full background see the
[PCM/ddX workflow page](https://open-quantum-platform.github.io/openqp-docs/workflows/pcm/)
and the [`[pcm]` keyword reference](https://open-quantum-platform.github.io/openqp-docs/keywords/pcm/).

## Input-file style

The runnable deck is [`inputs/h2o_pcm_ddpcm.oqp`](inputs/h2o_pcm_ddpcm.oqp). It
is an ordinary RHF single point with one extra call — `pcm(...)` — that switches
on the reaction field.

```text
rhf/6-31g* pcm                          # RHF reference in a ddPCM water continuum
geom="""
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
"""
```

Line by line:

- **`rhf/6-31g*`** is the reference wavefunction the reaction field couples to.
  Use `rhf` (closed shell) or `rohf` (open shell); **UHF is rejected**. There is no
  driver line, so the deck runs `energy()` — PCM/ddX is a single-point energy path,
  and `grad`/`opt` are outside its current scope.
- **`pcm`** turns the reaction field on. Written bare it uses the defaults —
  the ddX backend coupled to the reference SCF, the `ddpcm` model, and water
  (`epsilon=78.3553`) — which is the supported combination (energy driver,
  RHF/ROHF reference) the input checker enforces. Write options only to change
  them, for example `pcm(model=ddcosmo,epsilon=4.8)`; `ddlpb` is the third
  model, `solvent` a readable label, and `radii` the cavity-radii model
  (default `uff`).
- The default SCF threshold (`conv=1e-6`) is the right one for this energy-only
  coupling — the shipped examples note that a tighter 1e-8 gate is unreachable
  because the provisional coupling omits one small (~1e-7) Fock term, and 1e-6
  already reproduces the reference energy well within tolerance.

> **Note on `ispher`.** The reference deck in the OpenQP repo sets `ispher=true`
> to exercise spherical-harmonic AOs, but the docs state `ispher` is normally
> selected automatically from the basis convention for PCM/ddX inputs, so this
> minimal deck leaves it out.

## Python style

[`inputs/h2o_pcm_ddpcm.py`](inputs/h2o_pcm_ddpcm.py) produces the same solvated
energy through the OpenQP Python API. `job.theory.hf(...)` sets the reference and
`job.workflow.pcm(...)` maps one-to-one onto the `[pcm]` section (it also
pre-validates the ddX / reference_scf / RHF-or-ROHF scope before running).

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_pcm", silent=1)

job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)
job.theory.hf(basis="6-31g*")          # reference="rohf", multiplicity=2 for OH

job.workflow.pcm(
    enabled=True,
    backend="ddx",
    mode="reference_scf",
    model="ddpcm",
    epsilon=78.3553,
)

mol = job.run()
print("Solvated SCF energy (Hartree):", mol.get_scf_energy())
```

## Run it

Both styles run from the `inputs/` folder.

```bash
cd pcm-solvation/inputs

# Input-file style
openqp h2o_pcm_ddpcm.oqp

# Python style
python h2o_pcm_ddpcm.py
```

The `.oqp` run writes `h2o_pcm.log` (and a JSON summary); the Python script
prints the solvated SCF energy to stdout.

## Reading the output

The number to look at is the **SCF total energy** — it already contains the PCM
reaction-field contribution, so it is the *solvated* energy, not a gas-phase
one. In the log it is the converged SCF energy; from Python it is
`mol.get_scf_energy()` (or `mol.get_results()["energy"]`). The reference JSON
shipped with the repo example records `-76.02578997875146` Hartree for this
water case, so a converged run should land right there.

To see the actual **solvation effect**, run the same molecule twice — once with
`[pcm] enabled=false` (or with no `[pcm]` section / no `job.workflow.pcm(...)`
call) and once with it on — and subtract. The difference is the reaction-field
stabilization: a polar solute in water is lowered by a few kcal/mol to tens of
kcal/mol depending on its charge distribution. Watch the SCF also take **a few
more iterations** to converge with PCM on, since the reaction field is updated
as the density relaxes.

## References / manual

- PCM/ddX workflow — [open-quantum-platform.github.io/openqp-docs/workflows/pcm](https://open-quantum-platform.github.io/openqp-docs/workflows/pcm/)
- `[pcm]` keyword reference — [open-quantum-platform.github.io/openqp-docs/keywords/pcm](https://open-quantum-platform.github.io/openqp-docs/keywords/pcm/)
- `[scf]` keyword reference — [open-quantum-platform.github.io/openqp-docs/keywords/scf](https://open-quantum-platform.github.io/openqp-docs/keywords/scf/)
- Driving OpenQP from Python — [open-quantum-platform.github.io/openqp-docs/python-scripting](https://open-quantum-platform.github.io/openqp-docs/python-scripting/)
- PCM and domain-decomposition ddX literature — [references page](https://open-quantum-platform.github.io/openqp-docs/references/#pcm-and-ddx)
- OpenQP — [10.1021/acs.jctc.4c01117](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)
