# Quantum Krylov and quantum preconditioning

Martin Stoll and John W. Pearson.

| folder | content | status |
|---|---|---|
| `preconditioning/` | Paper "Preconditioning for quantum linear solvers in PDE-constrained optimization": `paper/` (LaTeX), `code/` (numerical checks, `qp_*.py`), `next_steps_preconditioning.md` (to-do list) | active |
| `krylov_note/` | Note "Quantum Krylov subspace methods" for numerical linear algebraists (`quantum_krylov_note.tex`), its experiments in `code/`, and `dead_ends.md` (ideas that did not work) | done |
| `archive/rational_krylov/` | Rational Krylov / Loewner idea, pilot and prior-art check | parked |

**Python.** Use the local virtual environment, `.venv/bin/python` (numpy, scipy, matplotlib; not committed). Scripts in a `code/` folder import each other, so keep them together.

**LaTeX.** Run `latexmk -pdf <file>.tex` inside the folder of the `.tex` file, then `latexmk -c`.
