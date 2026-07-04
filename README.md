# OpenQP Tutorials

Hands-on, runnable tutorials for **[OpenQP](https://github.com/Open-Quantum-Platform/openqp)**
(Open Quantum Platform). Every tutorial shows the calculation in **both** the
classic **input-file** style and the compact **Python API**.

📖 **Read the tutorials book:** <https://open-quantum-platform.github.io/openqp-tutorials/>

This repository is a [MkDocs](https://www.mkdocs.org/) site (same look and feel as
the [OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/)). Each
tutorial lives in `docs/<topic>/` with an `index.md` walkthrough and an `inputs/`
folder of ready-to-run decks.

## Install

```bash
pip install openqp      # the QM engine + `openqp` CLI
pip install openmm      # optional MM backend (needed for the QM/MM tutorials)
```

## Run a tutorial (either style)

```bash
openqp docs/mrsf-tddft/inputs/h2o_mrsf.inp     # input-file style
python docs/mrsf-tddft/inputs/h2o_mrsf.py      # Python-API style
```

## Contents

**Electronic structure** — Hartree-Fock & DFT · MP2 · TDDFT/TDHF · Spin-flip
TDDFT · MRSF-TDDFT · UMRSF-TDDFT
**Excited states & dynamics** — Spin-orbit coupling · Conical intersections ·
SOC-NAMD-QMMM
**Geometry & properties** — Optimization & TS · Hessians/frequencies/IR/Raman ·
NMR shielding · Population, moments & MRSF analysis
**Environment** — PCM/ddX solvation · ESPF QM/MM embedding
**SCF & basis** — SCF convergence & guesses · Effective core potentials

## Build the book locally

```bash
pip install -r docs/requirements.txt
mkdocs serve       # live preview at http://127.0.0.1:8000
```

## Contributing

A tutorial is `docs/<topic>/index.md` (motivation → theory → annotated input →
Python → run → interpret) plus an `inputs/` folder of small, fast decks. Add it to
the `nav` in `mkdocs.yml`. Link to the
[manual](https://open-quantum-platform.github.io/openqp-docs/) for reference depth
rather than duplicating it.

## License

Released under the same license as OpenQP (GPLv3).
