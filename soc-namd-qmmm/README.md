# SOC-NAMD-QMMM: excited-state dynamics with intersystem crossing in an environment

This tutorial walks you through a **spin-orbit-coupled nonadiabatic molecular
dynamics simulation with QM/MM embedding** — "SOC-NAMD-QMMM" — from the physics
behind the acronym to a deck you can run and the numbers it prints.

By the end you will understand:

- what each ingredient (**MRSF-TDDFT**, **NAMD**, **SOC**, **QM/MM**) contributes
  and why you would combine them,
- how the pieces map onto the `[scf]`, `[tdhf]`, `[md]`, and `[qmmm]` sections of
  an OpenQP input,
- how to run the calculation and read the trajectory it produces, and
- how to grow the minimal demo into a production run.

The runnable inputs are in [`inputs/`](inputs/); the reference manual for every
keyword is the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/),
in particular the [SOC-NAMD-QMMM workflow page](https://open-quantum-platform.github.io/openqp-docs/workflows/soc-namd-qmmm/)
and the [`[md]`](https://open-quantum-platform.github.io/openqp-docs/keywords/md/)
/ [`[qmmm]`](https://open-quantum-platform.github.io/openqp-docs/keywords/qmmm/)
keyword pages.

---

## 1. The problem, and why it needs four ingredients

Imagine a molecule that absorbs light, is promoted to an excited singlet state,
and then — instead of simply fluorescing back down — **crosses into a triplet
state** and does its chemistry from there. This *intersystem crossing* (ISC) is
at the heart of phosphorescence, photosensitizers, triplet-triplet annihilation,
and much of photobiology. Simulating it from first principles, *in a realistic
environment*, needs four things working together. Let's build up the acronym one
ingredient at a time.

### MRSF-TDDFT — the excited states

You need excited-state energies and, crucially, **gradients** (forces) to move
the nuclei. Ordinary TDDFT struggles near conical intersections and for
open-shell/diradical character. **MRSF-TDDFT** (Mixed-Reference Spin-Flip TDDFT)
fixes this: it spin-flips out of a **high-spin (triplet) reference** and mixes two
reference determinants, which removes the spin contamination of plain spin-flip
TDDFT and gives balanced singlet *and* triplet states — exactly what you need when
those states are about to cross. In the deck this is `[scf] multiplicity=3
type=rohf` (the triplet reference) plus `[tdhf] type=mrsf`.

### NAMD — letting the nuclei hop between surfaces

Once a molecule is excited, the Born-Oppenheimer approximation (nuclei on a
single fixed electronic surface) breaks down: the wavepacket can *switch*
electronic states as it moves. **Nonadiabatic molecular dynamics** models this.
OpenQP uses **Tully's fewest-switches surface hopping (FSSH)**: the nuclei move
classically on one "active" state, an electronic amplitude is propagated
alongside, and at each step the trajectory may *hop* to another state with a
probability set by the nonadiabatic coupling. `runtype=namd` selects this; the
`[md]` section controls the integrator, the active state, and the hopping.

### SOC — turning on singlet ↔ triplet transitions

Plain NAMD only hops between states of the **same** spin (internal conversion,
S→S). A singlet cannot become a triplet without the relativistic **spin-orbit
coupling** operator, which mixes spin. Turn it on with `[md] soc=True`. OpenQP
then propagates on a **spin-adiabatic ("SHARC-like") manifold** built from the
singlets *and* the individual triplet sublevels — for `nstate` singlets and `nt`
triplets that is `ns + 3·nt` states (each triplet contributes 3 M_s sublevels).
A hop between a singlet block and a triplet block **is** an intersystem crossing,
driven by the computed SOC matrix elements. This is what makes it *SOC*-NAMD.

### QM/MM — putting the molecule in its environment

Photochemistry rarely happens in vacuum. **QM/MM** treats the interesting
molecule quantum-mechanically (QM) and its surroundings — solvent, protein,
crystal — with a cheap classical force field (MM). OpenQP couples the two with
**ESPF electrostatic embedding**: the MM point charges polarize the QM density
through the electrostatic-potential-fitted (ESPF) operator, and the QM density
reacts back on the MM atoms, with an analytic, energy-conserving gradient.
`qmmm_flag=True` plus the `[qmmm]` section (a PDB, a force field, and which atoms
are QM) turn this on; the MM engine is [OpenMM](https://openmm.org).

**Put together**, SOC-NAMD-QMMM = excited-state surface-hopping dynamics of an
MRSF-TDDFT chromophore, *with* singlet↔triplet intersystem crossing, *embedded*
in an explicit MM environment.

---

## 2. The system for this tutorial

We use **formaldehyde (H₂CO) in a small cluster of 5 water molecules**. It is the
smallest system that shows every ingredient:

- Formaldehyde's low-lying **n→π\*** singlet (S₁) and triplet (T₁) states lie
  close in energy — a classic ISC playground.
- 5 waters give a real, polarizing MM environment without being slow.
- It runs in seconds, so you can iterate.

The QM region is the 4 formaldehyde atoms (`qm_atoms=0-3` in the PDB); the 5
waters are MM.

> **Note.** The shipped deck uses `[md] nstep=1` — a *single* nuclear step — so
> the example finishes instantly and can be used as a smoke test. To see actual
> dynamics and hops you raise `nstep` (see §6).

---

## 3. Prerequisites

```bash
pip install openqp        # the QM engine + CLI
pip install openmm        # the MM backend (required for QM/MM)
```

Check both import:

```bash
python -c "import oqp; import openmm; print('ok')"
```

If OpenMM is missing, QM/MM decks are reported **SKIPPED** rather than run.

---

## 4. The input, section by section

Open [`inputs/h2co-water_soc-namd-qmmm.inp`](inputs/h2co-water_soc-namd-qmmm.inp).
It has six sections. Here is what each one is doing.

### `[input]` — the QM subsystem and the run type

```ini
system=
   6   0.000000   0.000000   0.000000     # C
   8   0.000000   0.000000   1.203000     # O
   1   0.000000   0.943000  -0.589000     # H
   1   0.000000  -0.943000  -0.589000     # H
charge=0
runtype=namd            # nonadiabatic dynamics
basis=6-31g*
functional=bhhlyp
method=tdhf
qmmm_flag=True
```

`system` is the **QM geometry**, in the same order as `qmmm_flag`'s selection.
`runtype=namd` picks surface-hopping dynamics; `qmmm_flag=True` says "embed this
in an MM environment (configured in `[qmmm]`)". `bhhlyp` is a common half-and-half
functional for MRSF.

### `[scf]` — the high-spin reference

```ini
multiplicity=3
type=rohf
```

MRSF-TDDFT is built *on top of* a **triplet ROHF** reference. This is not the
state you care about — it is the mathematical starting point that MRSF spin-flips
from to reach the balanced singlet/triplet manifold.

### `[tdhf]` — the MRSF excited states

```ini
type=mrsf
nstate=2                # solve 2 MRSF singlet roots
multiplicity=3          # reference multiplicity (matches [scf])
```

`nstate=2` requests two singlet roots (S₀, S₁). With `soc=True` the code also
forms the triplets and their sublevels automatically, giving the spin-adiabatic
manifold described in §1.

### `[properties]`

```ini
grad=1                  # analytic gradients: needed to move the nuclei
```

### `[md]` — the dynamics and the hopping

This is the heart of a NAMD run.

```ini
nstep=1                 # number of nuclear steps
dt=0.25                 # nuclear time step (fs)
active=5                # initially populated spin-adiabatic state (0-based)
substep=50              # electronic sub-steps per nuclear step
init_temp=300           # sample initial velocities at 300 K (Maxwell)
velocity=maxwell
seed=3                  # fixed RNG seed -> reproducible hops
decoherence=edc         # energy-based decoherence (Granucci-Persico 2007)
trivial=True            # detect trivial (weakly-avoided) crossings
soc=True                # spin-orbit coupling ON -> allows ISC
thrshe=0.1              # gap gate: suppress spurious hops to S0
```

Key ideas:

- **`active`** is *which* state the trajectory starts on. On the spin-adiabatic
  manifold the states are ordered by energy across singlets and triplet
  sublevels, so `active=5` starts partway up the manifold (a bright singlet,
  here) rather than the ground state.
- **`substep`** matters because the electronic amplitude oscillates much faster
  than the nuclei move; the time-derivative couplings are integrated on this finer
  grid between nuclear steps.
- **`decoherence=edc`** corrects a well-known FSSH pathology (over-coherence) so
  populations relax physically. **`trivial=True`** stops the trajectory from
  missing a hop at a very sharp, weakly-avoided crossing.
- **`soc=True`** is the switch that turns *NAMD* into *SOC-NAMD*: without it you
  get internal conversion only; with it, singlet↔triplet ISC is possible.
- **`thrshe`** is a safety gate: near the Franck-Condon point the default would
  permit unphysical hops down to S₀; `0.1` Hartree is the recommended value.

### `[qmmm]` — the environment

```ini
pdb_file=formaldehyde_water.pdb          # full QM+MM system
forcefield_files=formaldehyde.xml tip3p.xml
qm_atoms=0-3                             # 0-based PDB indices that are QM
cutoff=NoCutoff                          # isolated cluster
embedding=electrostatic                  # full ESPF embedding
```

`pdb_file` holds *all* atoms (QM + MM); `qm_atoms` carves out the QM region;
everything else is MM, parameterized by the `forcefield_files`. `formaldehyde.xml`
only needs to supply the QM atoms' Lennard-Jones parameters — their electrostatics
come from ESPF, not from fixed MM charges. `cutoff=NoCutoff` is for an isolated
cluster; for a **solvated periodic box** you would use `cutoff=PME`, which turns
on the particle-mesh-Ewald branch of the embedding.

> The QM region here is a *whole molecule*. Cutting a covalent bond (carving a
> fragment out of a larger molecule) is a **covalent boundary** — supported by the
> single-point and ground-state QM/MM paths, but not by `runtype=namd`. See the
> [covalent-boundary docs](https://open-quantum-platform.github.io/openqp-docs/keywords/qmmm/#covalent-qmmm-boundaries).

---

## 5. Running it

From the `inputs/` folder (so the PDB/force-field files resolve):

```bash
cd soc-namd-qmmm/inputs
openqp h2co-water_soc-namd-qmmm.inp
```

What OpenQP does at **each nuclear step**:

1. Build the MM electrostatic potential at the QM atoms (via OpenMM) and hand it
   to the ESPF operator.
2. Run the embedded MRSF-TDDFT calculation → singlet and triplet energies and the
   **SOC** matrix elements between them.
3. Form the spin-adiabatic states and the gradient of the active state.
4. Integrate the electronic amplitude over `substep` sub-steps using the
   time-derivative couplings, and evaluate the fewest-switches hopping
   probabilities (now including S↔T channels because `soc=True`).
5. Possibly **hop** (with decoherence + trivial-crossing handling), then advance
   the nuclei by `dt` under the active-state force (QM atoms by the QM/ESPF force,
   MM atoms by the MM + coupling force).

The run writes a log (`<project>.log`) and a trajectory.

---

## 6. Reading the output and going further

The `.log` records, per step: the electronic-state energies, the **active state**,
the hopping probabilities, the SOC couplings, and the kinetic/potential/total
energy. A hop shows up as the active-state index changing between steps; an
**intersystem crossing** is a hop whose old and new states belong to different
spin blocks.

To turn the smoke test into a real simulation:

| Want to… | Change |
| --- | --- |
| See actual dynamics/hops | Raise `[md] nstep` (e.g. 500–2000) |
| Start on a different state | Change `[md] active` |
| Sample many trajectories | Vary `[md] seed` and average |
| Model a solvated box | `[qmmm] cutoff=PME` with a periodic PDB |
| More/other excited states | Raise `[tdhf] nstate` |
| Internal conversion only | `[md] soc=False` |

**On the SOC gradient.** With `soc=True` the default hops on the
**spin-adiabatic** manifold using a weighted-MCH diagonal gradient. If you need
the *exact* active-root gradient in the molecular-Coulomb-Hamiltonian (MCH) basis,
set `[md] soc_basis=mch`. See the workflow page for the trade-offs.

---

## 7. Recap

You combined four ideas into one simulation:

- **MRSF-TDDFT** gave balanced singlet/triplet excited states and gradients,
- **NAMD (FSSH)** let the nuclei hop between electronic surfaces,
- **SOC** opened the singlet↔triplet (intersystem-crossing) channel, and
- **QM/MM (ESPF)** embedded the chromophore in an explicit environment.

That is SOC-NAMD-QMMM. From here, the [manual](https://open-quantum-platform.github.io/openqp-docs/workflows/soc-namd-qmmm/)
gives the full input contract, the periodic (PME) setup, and the compact
`job.qmmm(...)` / `job.workflow.namd(...)` Python API.

## References

- Mixed-Reference Spin-Flip TDDFT — [10.1021/acs.jpclett.3c02296](https://doi.org/10.1021/acs.jpclett.3c02296)
- Relativistic (SOC) MRSF-TDDFT — [10.1021/acs.jctc.2c01036](https://doi.org/10.1021/acs.jctc.2c01036)
- ESPF QM/MM embedding — Huix-Rotllant & Ferré, [10.1021/acs.jctc.0c01075](https://doi.org/10.1021/acs.jctc.0c01075)
- Energy-based decoherence — Granucci & Persico, J. Chem. Phys. 126, 134114 (2007)
- OpenQP — [10.1021/acs.jctc.4c01117](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)
