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

The runnable deck is [`inputs/h2o_pcm_ddpcm.inp`](inputs/h2o_pcm_ddpcm.inp). It
is an ordinary RHF single-point with one extra section — `[pcm]` — that switches
on the reaction field.

```ini
# Water RHF energy in a ddPCM dielectric continuum (energy-only PCM/ddX).
[input]
system=
   8   0.000000000   0.000000000  -0.041061554   # O
   1  -0.533194329   0.533194329  -0.614469223   # H
   1   0.533194329  -0.533194329  -0.614469223   # H
charge=0
runtype=energy        # PCM/ddX is an energy-only path
basis=6-31g*
method=hf             # HF reference (Kohn-Sham DFT also allowed)

[guess]
type=huckel           # extended-Huckel initial orbitals
save_mol=false

[scf]
multiplicity=1        # closed shell
type=rhf              # RHF (or rohf for open-shell); UHF is not supported
conv=1.0e-6           # 1e-6 is the recommended SCF gate for this path

[pcm]
enabled=true          # turn the reaction field on
backend=ddx           # domain-decomposition solver (only ddx is implemented)
mode=reference_scf    # couple PCM to the reference SCF energy
model=ddpcm           # ddcosmo | ddpcm | ddlpb  (ddx continuum flavors)
epsilon=78.3553       # static dielectric constant (water); must be > 1
```

Section by section:

- **`[input]`** — geometry, charge, basis, and `method=hf`. Crucially
  `runtype=energy`: PCM/ddX is a single-point energy path, so `grad`/`optimize`
  runtypes are outside its current scope.
- **`[guess]`** — a Huckel starting guess; `save_mol=false` keeps the demo from
  leaving a restart file.
- **`[scf]`** — the reference wavefunction the reaction field couples to. Use
  `type=rhf` (closed shell) or `type=rohf` (open shell); **UHF is rejected**.
  `conv=1.0e-6` is the deliberate SCF threshold for this energy-only coupling —
  the shipped examples note that a tighter 1e-8 gate is unreachable because the
  provisional coupling omits one small (~1e-7) Fock term, and 1e-6 already
  reproduces the reference energy well within tolerance.
- **`[pcm]`** — the solvation controls. `enabled=true` activates it;
  `backend=ddx` + `mode=reference_scf` + `runtype=energy` + an RHF/ROHF
  reference is the supported combination the input checker enforces. `model`
  picks the continuum flavor (swap in `ddcosmo` or `ddlpb`), and `epsilon` sets
  the solvent's dielectric constant. Optional keys not shown here: `solvent`
  (a readable label), `radii` (cavity-radii model, default `uff`).

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
openqp h2o_pcm_ddpcm.inp

# Python style
python h2o_pcm_ddpcm.py
```

The `.inp` run writes `h2o_pcm.log` (and a JSON summary); the Python script
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
