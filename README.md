# OpenQP Tutorials

Hands-on, runnable tutorials for **[OpenQP](https://github.com/Open-Quantum-Platform/openqp)**
(Open Quantum Platform) — the MRSF-TDDFT quantum-chemistry package.

Where the pieces fit:

| Resource | What it is |
| --- | --- |
| **[openqp](https://github.com/Open-Quantum-Platform/openqp)** | the code: the QM engine, the `openqp` CLI, and the packaged `examples/`. |
| **[openqp-docs](https://open-quantum-platform.github.io/openqp-docs/)** | the manual: reference documentation for every method, workflow, and keyword. |
| **openqp-tutorials** (this repo) | guided, explanatory walkthroughs — the *why* and *how*, end to end, with decks you can run. |

Tutorials teach a workflow from motivation to result; the manual is the reference
you reach for once you know what you are doing.

## Tutorials

| Tutorial | What you learn |
| --- | --- |
| [SOC-NAMD-QMMM](soc-namd-qmmm/) | Excited-state surface-hopping dynamics of an MRSF-TDDFT chromophore, with spin-orbit **intersystem crossing**, embedded in an explicit MM environment — built up one ingredient at a time and run on formaldehyde in water. |

*(more coming — contributions welcome)*

## Getting started

```bash
pip install openqp      # the QM engine + `openqp` CLI
pip install openmm      # optional MM backend (needed for the QM/MM tutorials)
```

Each tutorial folder has a `README.md` walkthrough and an `inputs/` directory
with ready-to-run decks:

```bash
cd soc-namd-qmmm/inputs
openqp h2co-water_soc-namd-qmmm.inp
```

## Contributing

A tutorial is a self-contained folder with a `README.md` (motivation → theory →
annotated input → run → interpret → extend) and an `inputs/` directory of small,
fast, runnable decks. Prefer minimal systems that finish in seconds so readers can
iterate. Link to the [manual](https://open-quantum-platform.github.io/openqp-docs/)
for reference depth rather than duplicating it.

## License

Released under the same license as OpenQP (GPLv3).
