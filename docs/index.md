# OpenQP Tutorials

Hands-on, runnable tutorials for **[OpenQP](https://github.com/Open-Quantum-Platform/openqp)**
(Open Quantum Platform) — the MRSF-TDDFT quantum-chemistry package. Every tutorial
shows the calculation in **both styles**: the concise **`.oqp` input deck** and the
compact **Python API**.

Where the pieces fit:

| Resource | What it is |
| --- | --- |
| **[openqp](https://github.com/Open-Quantum-Platform/openqp)** | the code: the QM engine, the `openqp` CLI, and the packaged `examples/`. |
| **[openqp-docs](https://open-quantum-platform.github.io/openqp-docs/)** | the manual: reference documentation for every method, workflow, and keyword. |
| **openqp-tutorials** (this book) | guided, explanatory walkthroughs — the *why* and *how*, end to end. |

Tutorials teach a workflow from motivation to result; the manual is the reference
you reach for once you know what you are doing.

## Install

```bash
pip install openqp      # the QM engine + `openqp` CLI
pip install openmm      # optional MM backend (needed for the QM/MM tutorials)
```

## How to use

Each tutorial page walks through the physics, an **annotated `.oqp` deck**, the
**equivalent Python script**, how to run both, and how to read the output. The
runnable files live next to each tutorial in its `inputs/` folder. Every deck
uses a small, fast system (water, ethylene, formaldehyde) so you can iterate in
seconds.

Run a tutorial either way:

```bash
openqp <tutorial>/inputs/<deck>.oqp     # input-file style
python <tutorial>/inputs/<deck>.py      # Python-API style
```

## The `.oqp` input format

Every deck in this book is written in OpenQP's concise `.oqp` format. A deck is
built from four kinds of item — the route, the driver, and any options on one
line, followed by the geometry — and it says only what differs from the
defaults:

```text
mrsf(nstate=3)/bhhlyp/6-31g* grad(S1)
geom="h2o.xyz"
```

| Item | Example | What it does |
| --- | --- | --- |
| **route** (always first) | `mrsf(nstate=3)/bhhlyp/6-31g*` | names the physical model, the functional, and the basis. Model options go in parentheses. |
| **driver** (at most one) | `grad(S1)` | names the calculation and its target state. `energy()` is the default. |
| **options / section calls** | `charge=1`, `scf(conv=1e-10)` | top-level physical settings, plus exact legacy-section calls for anything the defaults do not cover. |
| **geometry** | `geom="h2o.xyz"` or an inline `geom` block | the molecule: an `.xyz`/`.pdb` path, or coordinates in triple quotes. |

Two things the format does deliberately:

- **States are physical.** You write `grad(S1)` or `meci(S0,S1)`; you never work
  out which internal response root that is. (Spin-flip roots are not spin-adapted
  before diagonalization, so SF decks use `root=N` instead.)
- **References are implied by the model.** `mrsf` means a high-spin triplet ROHF
  reference, `umrsf` a UHF one; you do not restate them.

Anything the concise surface does not name is still reachable through an exact
section call — `tdhf(nvdav=30)`, `dftgrid(rad_npts=96,ang_npts=302)` — so no
keyword is lost. The older sectioned `.inp` format is still read by `openqp`;
see the [manual](https://open-quantum-platform.github.io/openqp-docs/) for it.

## The tutorials

- **Electronic structure** — [Hartree-Fock and DFT](hf-and-dft/), [MP2](mp2/),
  [TDDFT and TDHF](tddft-and-tdhf/), [Spin-flip TDDFT](sf-tddft/),
  [MRSF-TDDFT](mrsf-tddft/), [UMRSF-TDDFT](umrsf-tddft/).
- **Excited states and dynamics** — [Spin-orbit coupling](spin-orbit-coupling/),
  [Conical intersections](conical-intersections/), [SOC-NAMD-QMMM](soc-namd-qmmm/).
- **Geometry and properties** — [Geometry optimization and TS](geometry-optimization/),
  [Hessians, frequencies, IR/Raman](vibrational-analysis/), [NMR shielding](nmr-shielding/),
  [Population, moments, MRSF analysis](properties-and-population/).
- **Environment** — [PCM/ddX solvation](pcm-solvation/),
  [ESPF QM/MM embedding](qmmm-embedding/).
- **SCF and basis** — [SCF convergence and guesses](scf-convergence/),
  [Effective core potentials](effective-core-potentials/).

New tutorials welcome — see the [repository](https://github.com/Open-Quantum-Platform/openqp-tutorials).
