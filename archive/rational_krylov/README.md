# Archive: rational Krylov / Loewner methods for quantum eigenvalue problems

Parked on 2026-11-03. Kept for reference only; nothing here is maintained.

**Idea.** Resolvent filters (z_k − H)⁻¹v₀ give a Loewner pencil built from r(z) = v₀*(z − H)⁻¹v₀. r(z) is estimated as a Laplace transform of quantum time-series data, by random-time Monte Carlo.

**Why parked.** The core mechanism is classical filter diagonalization, and it has already been brought to quantum computers. The closest work is quantum prolate diagonalization (QPD). See `prior_art.md`. What might still be new is narrow:
- a continuous-time Laplace/Loewner formulation;
- a shared-shot estimator;
- NLA theory;
- decoherence treated as a pole shift.

A head-to-head comparison with QPD would be required.

**Files**
- `manuscript_outline.md`: planned paper, claimed advantages, pilot results, risks.
- `prior_art.md`: prior-art check (filter diagonalization, QPD, quantum inverse iteration, fault-tolerant contour methods).
- `qk_rational_pilot.py`, `fig_rational_pilot.pdf`: pilot experiments E1–E3 (noiseless, equal shots, damping).
- `qk_rational_fair.py`: **unfinished**. A "fair" comparison (damping plus shot noise, equal total evolution time, two-stage pole choice) written by an agent that was stopped before it produced results. Never validated.
