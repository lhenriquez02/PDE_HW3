# Numerical Solutions to PDEs — Working with Immersed Boundaries

Mini-project implementing the **Immersed Boundary Method (IBM)** — technically
the **Force Coupling Method (FCM)**, which differs only in the choice of kernel
(a 2D Gaussian regularized delta with width a = 1.2Δx). The 2D Poisson
finite-difference code is extended to handle an *internal* Dirichlet condition
by introducing spreading/interpolation operators and solving the resulting
saddle-point system with **GMRES**.

A written summary of the whole assignment, including every figure referenced
below, lives in `report.tex` / `report.pdf`.

## Problems

### Problem 1 — Poisson Equation (sanity check)
Solve ∇²u = f on [0, 2π]² with Dirichlet conditions on all sides, using the
five-point finite-difference Laplacian, for the manufactured solution
u(x, y) = cos(2x) sin(2y). Second-order accuracy is verified: the error table
and the log–log error plot show the max-norm error proportional to Δx².

*Code:* `Code/question5_hw2_problem1.m` · *Figures:* `Graphics/poisson_hw2_problem1_*`

### Problem 2 — Poisson Equation with an Immersed Boundary
Same setup, but with an additional internal Dirichlet condition imposed on a
circular immersed boundary (center (π, π), radius 1):

- x_ib = π + cos(θ), y_ib = π + sin(θ), θ ∈ [0, 2π]
- u(x_ib, y_ib) = cos(2x_ib) sin(2y_ib)

The test solution is unchanged, so it can still be used for verification.
Procedure:

1. **Discretize** the IB curve with N_ib points spaced roughly equal to Δx of
   the FD grid.
2. **Spread:** place a charge q_k on each IB point. The charges modify the PDE
   to ∇²u = f + Σₖ qₖ δ_a(x − x_k^ib, y − y_k^ib) = f + S^ib(q), where δ_a is a
   2D Gaussian with variance a = 1.2Δx.
3. **Interpolate:** enforce the internal boundary condition by interpolating u
   onto the IB grid, J^ib(u) = cos(2x_ib) sin(2y_ib), giving enough equations
   for the unknowns.
4. **Solve** the saddle-point system

   ```
   [ ∇²   −S^ib ] [ u ]   [ f ]
   [ J^ib    0  ] [ q ] = [ cos(2x_ib) sin(2y_ib) ]
   ```

   matrix-free with **GMRES** (no explicit matrices; block preconditioner built
   from `L` and the Schur complement `J L⁻¹ S`).
5. **Verify convergence:** the error decreases at every refinement and, in the
   discrete L² norm, approaches the promised **O(Δx^1.5)** reduced order. The
   max-norm error instead plateaus: the near-boundary charge layer of the
   regularized delta does not vanish in the pointwise norm, and the GMRES
   tolerance doubles as a regularization knob — solving the saddle system "too
   exactly" amplifies near-null charge modes and hurts the pointwise error.
6. The result looks identical to the original solution (see the report
   figures).

*Code:* `Code/problem2.m` · *Figures:* `Graphics/poisson_hw2_problem2_*`

### Problem 3 — An Almost Minimal Surface
Same immersed boundary setup with f = 0 and homogeneous Dirichlet conditions on
the box: ∇²u = S^ib(q), u = 0 on ∂Ω, u(x_ib, y_ib) = cos(2x_ib) sin(2y_ib).
There is no manufactured exact solution; the resulting shape is closely related
to a *minimal surface* — the form a soap film would take on the given
wire-frame.

*Code:* `Code/problem3.m` · *Figures:* `Graphics/poisson_hw2_problem3_*`

### Extension — Multiple Shapes (two circles)
The immersed-boundary machinery generalizes to several disjoint curves at once.
Two circles of radius 0.65 centered at (π ± 1.25, π) are held at u = +1 and
u = −1; the marker lists are merged into a single saddle-point system and solved
with the same GMRES solver. The left interior fills to +1, the right to −1, with
a dipole-like decay through the gap and out to the zero box boundary.

*Code:* `Code/two_circles.m` · *Figures:* `Graphics/poisson_hw2_two_circles_*`

## Plotting conventions
All figures follow the class textbook's visual-communication chapter:

- the colorblind-safe class palette (teal/magenta family), with a divergent
  colormap (teal → charcoal → magenta) for oscillatory fields,
- discrete data plotted with symbols; reference slopes drawn as black dashed
  lines, kept out of the legend and labeled directly on the plot,
- LaTeX (Latin Modern) labels, ticks and legends at large font sizes,
- log–log axes for error/convergence plots.

## Repository layout

- `report.tex`, `report.pdf` — written report, including all figures
- `Code/` — MATLAB sources: `question5_hw2_problem1.m`, `problem2.m`,
  `problem3.m`, `two_circles.m`
- `Graphics/` — the exported figure PDFs referenced by the report
- `PDES_IB_HW2.pdf` — assignment handout (the only gitignored file)

Build the report with `latexmk report.tex`. The figure PDFs in `Graphics/`
are produced by the MATLAB scripts (`exportgraphics`).
