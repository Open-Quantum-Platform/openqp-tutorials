# OpenQP Tutorials

Hands-on, runnable tutorials for **[OpenQP](https://github.com/Open-Quantum-Platform/openqp)**
(Open Quantum Platform) — the MRSF-TDDFT quantum-chemistry package. Every tutorial
shows the calculation in **both styles**: the classic **input-file** deck and the
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

Each tutorial page walks through the physics, an **annotated input deck**, the
**equivalent Python script**, how to run both, and how to read the output. The
runnable files live next to each tutorial in its `inputs/` folder. Every deck
uses a small, fast system (water, ethylene, formaldehyde) so you can iterate in
seconds.

Run a tutorial either way:

```bash
openqp <tutorial>/inputs/<deck>.inp     # input-file style
python <tutorial>/inputs/<deck>.py      # Python-API style
```

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
