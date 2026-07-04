# Population, moments, and MRSF excited-state analysis

A single-point energy tells you *how stable* a molecule is; the **population and
moment analysis** on top of it tells you *where the electrons went*. This
tutorial adds four opt-in property analyses to an ordinary RHF calculation of
water — the **electric dipole moment** and three flavours of **atomic partial
charge** (Mulliken, Lowdin, and RESP) — and shows how to read each number both
from the log and from the Python result dictionary. These are the same
descriptors you reach for when you need to assign a formal charge distribution,
parameterize a force field, or sanity-check a dipole against experiment. The
machinery is entirely opt-in: nothing beyond the energy is computed unless you
list it in `[properties] scf_prop`.

## A little theory

An SCF calculation converges a one-electron density. Every property in this
tutorial is a different *contraction* of that density with a set of operators:

- **Electric dipole (`el_mom`)** — the first electric moment, the expectation
  value of the position operator weighted by the charge density plus the nuclear
  contribution. It is the leading term describing how the molecule couples to a
  uniform external field, and it vanishes only for centrosymmetric systems.
- **Mulliken charges (`mulliken`)** — partition the density over atoms by
  splitting each overlap population *equally* between the two atoms of a basis-
  function pair. Cheap and intuitive, but famously basis-set sensitive.
- **Lowdin charges (`lowdin`)** — the same idea after a symmetric (S^-1/2)
  orthogonalization of the basis, which removes the arbitrary off-diagonal
  overlap term and usually gives more stable numbers than Mulliken.
- **RESP charges (`resp`)** — **R**estrained **E**lectro**S**tatic **P**otential
  charges: fit atom-centered point charges so that they reproduce the *quantum*
  electrostatic potential on a grid of points around the molecule. These are the
  charges you want for classical force fields, because they reproduce the field a
  neighbour actually feels rather than an orbital-partitioning convention.

Because Mulliken and Lowdin are basis-set–dependent partitions of the same
density while RESP is fit to a physical observable (the ESP), the three will not
agree exactly — comparing them is itself a useful diagnostic. For the operator
definitions and the population-analysis background, see the
[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/).

## Input-file style

The runnable deck is [`inputs/h2o_properties.inp`](inputs/h2o_properties.inp) —
water in the 6-31G\* basis, a closed-shell RHF reference, with all four property
analyses switched on. Annotated:

```ini
[input]
system=
   8   0.000000000   0.000000000  -0.041061554   # O   (Angstrom)
   1  -0.533194329   0.533194329  -0.614469223   # H
   1   0.533194329  -0.533194329  -0.614469223   # H
charge=0
runtype=energy          # single-point; properties ride on top of the energy
basis=6-31g*
method=hf               # Hartree-Fock reference (no functional -> pure HF)

[guess]
type=huckel             # extended-Huckel initial orbital guess

[scf]
multiplicity=1          # singlet ground state
type=rhf                # closed-shell restricted HF reference

[properties]
scf_prop=el_mom,mulliken,lowdin,resp   # the opt-in analyses to run
```

Key points:

- **`[properties] scf_prop` is the whole tutorial.** It is a comma-separated list;
  each token triggers one analysis after the SCF converges. The accepted values
  are `el_mom`, `mulliken`, `lowdin`, `resp`, and `nmr` (NMR shielding is covered
  in its own tutorial). Order does not matter, and you can request any subset.
- **Only what you list is surfaced.** `scf_prop` defaults to empty, so a bare
  `[input]`/`[scf]` deck prints just the energy. Requesting a property is what
  makes its numbers appear in the log and in the JSON/Python results — this is
  deliberate, so the regression-tested output contains exactly the descriptors
  you asked for.
- **The reference is an ordinary SCF**, chosen in `[scf] type` (`rhf` here). Water
  is closed-shell so `rhf` is natural; the analyses work the same on a `uhf` or
  `rohf` reference for open-shell systems — you would just change `type` and
  `multiplicity`. Because `method=hf` and no `functional` is given, this is pure
  Hartree-Fock; set a `functional` to run the identical property analysis on a
  DFT density instead.
- **`runtype=energy`** — these are ground-state SCF properties evaluated at the
  input geometry; no gradient is needed.

## Python style

The equivalent calculation with the OpenQP Python API is
[`inputs/h2o_properties.py`](inputs/h2o_properties.py). Two helper calls set up
the run: `job.theory.hf(...)` writes `[input] method=hf` and `[scf] type=rhf`,
and `job.settings.properties(...)` writes the raw `[properties]` section verbatim
— `job.settings.<section>(...)` is the general escape hatch for any OpenQP
section keyword that does not have a dedicated helper.

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_properties", silent=1)

# Same geometry (Angstrom) as the input deck.
job.molecule(
    """
O   0.000000000   0.000000000  -0.041061554
H  -0.533194329   0.533194329  -0.614469223
H   0.533194329  -0.533194329  -0.614469223
""",
    charge=0,
    multiplicity=1,
)

# RHF reference at 6-31G*. theory.hf() sets [input] method=hf and [scf] type=rhf.
job.theory.hf(basis="6-31g*", reference="rhf")

# Opt-in SCF properties: [properties] scf_prop=el_mom,mulliken,lowdin,resp.
job.settings.properties(scf_prop="el_mom,mulliken,lowdin,resp")

mol = job.run()
results = mol.get_results()

print("SCF energy (Ha):", results["energy"])
print("Dipole (a.u.):", results["dipole"])
print("Mulliken charges:", results["mulliken_charges"])
print("Lowdin charges:", results["lowdin_charges"])
print("RESP charges:", results["resp_charges"])
```

The five `print` lines pull each requested property straight out of the result
dictionary. The keys mirror the `scf_prop` tokens: `el_mom` -> `dipole`,
`mulliken` -> `mulliken_charges`, `lowdin` -> `lowdin_charges`, `resp` ->
`resp_charges`. A key only appears if you requested it — drop a token from
`scf_prop` and the matching key disappears from `results`.

## Run it

Input-file style (from the `inputs/` folder):

```bash
cd properties-and-population/inputs
openqp h2o_properties.inp
```

Python style:

```bash
cd properties-and-population/inputs
python h2o_properties.py
```

Both need OpenQP installed (`pip install openqp`) and produce the same numbers.

## Reading the output

The SCF **energy** is the primary number; the four properties are the payload of
this tutorial. For this water / 6-31G\* / RHF example the validated run gives:

| Quantity | Value |
| --- | --- |
| `energy` (SCF total) | `-76.0107465151` Ha |
| `dipole` (z-component) | `-0.8651` a.u. (x, y ≈ 0 by symmetry) |
| `mulliken_charges` (O, H, H) | `-0.8688`, `+0.4344`, `+0.4344` |
| `lowdin_charges` (O, H, H) | `-0.7312`, `+0.3656`, `+0.3656` |
| `resp_charges` (O, H, H) | `-0.8746`, `+0.4361`, `+0.4385` |

- In the **log file** (`<project>.log`) look for the converged SCF energy, then
  the dipole vector and each atom-charge table printed after it.
- From **Python**, `results["energy"]` is the SCF total, `results["dipole"]` is
  the dipole as a 3-vector (a.u.), and `results["mulliken_charges"]`,
  `results["lowdin_charges"]`, `results["resp_charges"]` are per-atom lists in the
  same order as the input geometry (O, H, H). These are exactly the fields written
  to `<project>.json`.
- **Read them against each other.** All three charge sets agree on the sign
  (oxygen negative, hydrogens positive) and roughly on magnitude, but Mulliken and
  RESP put more charge on oxygen than Lowdin does — the expected spread between an
  orbital-partitioning charge and an ESP-fitted one. The dipole points along the
  C2 axis (the molecular z), with the in-plane components zero to numerical noise,
  which is the symmetry check you would want to see for water.

## Manual

- SCF properties and `[properties]` keyword reference (`scf_prop` values,
  population and moment analysis):
  <https://open-quantum-platform.github.io/openqp-docs/keywords/properties/>
- `[scf]` keyword reference (`type`, `multiplicity` — choosing the SCF
  reference): <https://open-quantum-platform.github.io/openqp-docs/keywords/scf/>
- Running OpenQP from Python (the `job.theory.hf(...)` and
  `job.settings.<section>(...)` idioms):
  <https://open-quantum-platform.github.io/openqp-docs/python-scripting/>
