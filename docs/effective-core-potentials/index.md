# Effective core potentials (ECP): replacing the core with a pseudopotential

Heavy atoms have a lot of core electrons that barely change during chemistry but
are expensive to describe and — past the third row — genuinely **relativistic**.
An **effective core potential (ECP)**, also called a pseudopotential, throws those
inert core electrons out of the explicit calculation and replaces them with an
analytic potential acting on the remaining valence electrons. You reach for one
whenever a molecule contains heavier elements (roughly Br and below, and anything
in the 4th row or heavier): the ECP shrinks the number of electrons and basis
functions *and* folds scalar-relativistic effects into the potential, so a
nonrelativistic code gives near-relativistic valence energetics. This tutorial
computes a BHHLYP energy of hydrogen bromide (HBr) with an ECP on bromine and an
ordinary all-electron basis on hydrogen.

## A little theory

An ECP splits the electrons of an atom into a frozen **core** and an active
**valence**. The core is not solved for; instead its effect on the valence
electrons — nuclear attraction screened by the core, plus the Pauli repulsion
that keeps valence orbitals orthogonal to the (absent) core — is baked into a
short, angular-momentum-projected radial potential:

```text
V_ECP(r) = V_local(r) + Σ_l V_l(r) |l><l|
```

Each `V_l` is a small sum of Gaussian terms fitted to reproduce
all-electron (often relativistic) reference calculations. Because the potential
carries the relativistic contraction of the core, an ECP is the standard cheap
route to scalar-relativistic valence chemistry for heavy atoms. In practice you
never write the potential yourself: **correlation-consistent "-PP" basis sets**
(e.g. `aug-cc-pVDZ-PP`) ship a matched valence basis *and* the ECP that goes with
it. OpenQP reads that ECP from the basis-set data and evaluates its integrals
through the [libecpint](https://github.com/robashaw/libecpint) library; the rest
of the SCF is unchanged. The one thing to get right is the **per-atom basis**:
heavy atoms get the `-PP` set (with its ECP), light atoms keep an all-electron
set. For depth, see the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/hbr_ecp_energy.inp`](inputs/hbr_ecp_energy.inp) —
HBr, a BHHLYP single point, with `aug-cc-pVDZ-PP` (ECP) on Br and the
all-electron `aug-cc-pVDZ` on H. Annotated:

```ini
[input]
system=
   35   0.000000000   0.000000000   0.000000000     # Br  (Z=35, Angstrom)
   1    0.000000000   0.000000000   1.407611        # H
charge=0
runtype=energy          # single-point energy
functional=bhhlyp       # half-and-half hybrid; DFT XC on top of the HF path
basis=aug-cc-pVDZ-PP;aug-cc-pVDZ   # per-atom basis, in atom order (see below)
method=hf               # the SCF engine (KS-DFT is driven through it)

[guess]
type=huckel             # extended-Huckel initial orbital guess
save_mol=false

[scf]
multiplicity=1          # closed-shell singlet
type=rhf                # restricted HF/KS reference

[dftgrid]
rad_type=becke          # Becke radial grid for the DFT XC quadrature
```

Key points:

- **`basis=aug-cc-pVDZ-PP;aug-cc-pVDZ`** is the whole trick. The semicolon-separated
  list assigns one basis **per atom, in the same order as `system`**: the first
  atom (Br) gets `aug-cc-pVDZ-PP`, the second (H) gets `aug-cc-pVDZ`. The `-PP`
  suffix ("pseudopotential") is what makes OpenQP load bromine's ECP and treat its
  core implicitly; hydrogen is light, so it stays all-electron. There is no
  separate "turn on ECP" keyword — **choosing a `-PP` basis for an atom is how you
  request its ECP**, and the potential is built automatically through libecpint.
- **`functional=bhhlyp`** requests the BHHLYP hybrid, so this is a DFT single
  point. `method=hf` names the underlying SCF engine that the Kohn-Sham build runs
  through; for a bare Hartree-Fock energy you would leave `functional` empty. The
  ECP mechanics are identical either way.
- **`[scf] type=rhf` with `multiplicity=1`** — HBr is a closed-shell singlet. Note
  the electron count is over the *valence*: with the Br core removed by the ECP,
  the SCF handles far fewer electrons than an all-electron HBr calculation would.
- **`[dftgrid] rad_type=becke`** only affects the numerical XC quadrature for the
  functional; it is unrelated to the ECP.

To run the same molecule **all-electron** (no ECP), you would drop the `-PP` and
give bromine an all-electron relativistic-recontracted set — but that reintroduces
all 35 electrons of Br and its tight core basis functions, which is exactly the
cost the ECP is there to avoid.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/hbr_ecp_energy.py`](inputs/hbr_ecp_energy.py). The ECP is requested the
same way — a **per-atom list** passed to `job.settings.basis(...)`, in atom order —
so the `-PP` entry lands on Br and the all-electron entry on H.

```python
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
```

The moving parts map one-to-one onto the input file:

- **`job.theory.dft(functional="bhhlyp")`** is the `[input] functional=bhhlyp`
  line (and the `method=hf` engine underneath) — a closed-shell BHHLYP single
  point.
- **`job.settings.basis([...])`** is the `basis=...;...` line. Passing a **Python
  list** is the per-atom form: `["aug-cc-pVDZ-PP", "aug-cc-pVDZ"]` matches the two
  atoms of `job.molecule(...)` in order. A single string (e.g.
  `job.settings.basis("aug-cc-pVDZ-PP")`) would instead apply one basis to every
  atom. As in the input file, the ECP comes along for free with the `-PP` entry —
  no extra call is needed.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd effective-core-potentials/inputs
openqp hbr_ecp_energy.inp
```

Python style:

```bash
cd effective-core-potentials/inputs
python hbr_ecp_energy.py
```

Both need OpenQP installed (`pip install openqp`) and compute the identical
energy. OpenQP builds the ECP integrals through libecpint, which is bundled with
the standard OpenQP build — no extra install step.

## Reading the output

This is an energy-only run, so the number to read is the **converged SCF (BHHLYP)
total energy**:

- From **Python**, `mol.get_scf_energy()` returns the total energy, and
  `mol.get_results()` gives the full results dictionary (the same fields written to
  `<project>.json`).
- In the **log file** (`hbr_ecp_energy.log`) look for the final converged SCF
  energy. A quick sanity check that the ECP is active is the **electron count**:
  with bromine's core removed the SCF reports far fewer electrons than the 36 a
  fully all-electron HBr would have — the missing electrons are exactly Br's frozen
  core. If instead you see the full all-electron count, the `-PP` basis was not
  picked up (usually a per-atom ordering slip in `basis=`).

Because the ECP folds scalar-relativistic core effects into the valence, the HBr
bond energetics from this cheap valence-only run are close to a full relativistic
all-electron treatment — at a fraction of the cost.

## References / manual

- OpenQP manual (basis sets, ECP / `-PP` usage, and the full keyword contract):
  <https://open-quantum-platform.github.io/openqp-docs/>
- Running OpenQP from Python (the `job.theory.dft(...)` and
  `job.settings.basis([...])` idioms):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
- libecpint, the ECP-integral library OpenQP uses:
  <https://github.com/robashaw/libecpint>
