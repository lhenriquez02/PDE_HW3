# Numerical Solutions to PDEs — Working with Immersed Boundaries

Mini-project implementing the **Immersed Boundary Method (IBM)** — technically the
**Force Coupling Method (FCM)**, which differs only in the choice of kernel. The 2D
Poisson finite-difference code from HW0 is extended to handle an *internal* Dirichlet
boundary condition by introducing spreading/interpolation operators and solving the
resulting saddle-point system with a **matrix-free GMRES** solver.

## Problems

### Problem 1 — Poisson Equation (sanity check)
Solve ∇²u = f(x, y) on [0, 2π] × [0, 2π] with Dirichlet boundary conditions on all
sides, using finite differences. Choose f and boundary data so that

> u(x, y) = cos(2x) sin(2y)

is the exact test solution. Verify the code converges with **second order** in Δx.

### Problem 2 — Poisson Equation with an Immersed Boundary
Same setup, but now with an additional internal Dirichlet condition imposed on a
circular immersed boundary:

- x_ib = π + cos(θ), y_ib = π + sin(θ), θ ∈ [0, 2π]
- u(x_ib, y_ib) = cos(2x_ib) sin(2y_ib)

The test solution is unchanged, so it can still be used for verification. Procedure:

1. **Discretize** the IB curve with N_ib points spaced roughly equal to Δx of the FD grid.
   - Experiment with the spacing: effects of IB points much farther apart / much closer
     together than Δx, on both accuracy and speed.
2. **Spread:** place a charge q_k on each IB point. The charges modify the PDE to
   ∇²u = f + Σₖ qₖ δ_a(x − x_k^ib, y − y_k^ib) = f + S^ib(q), where δ_a is a 2D Gaussian
   with variance a = 1.2Δx.
3. **Interpolate:** enforce the internal boundary condition by interpolating u onto the
   IB grid, J^ib(u) = cos(2x_ib) sin(2y_ib), giving enough equations for the unknowns.
4. **Solve the saddle-point system**
   ```
   [ ∇²   −S^ib ] [ u ]   [ f ]
   [ J^ib    0  ] [ q ] = [ cos(2x_ib) sin(2y_ib) ]
   ```
   - (a) Do **not** construct the matrix explicitly — write a matrix-free apply routine
     using the S^ib and J^ib functions, then solve with **GMRES**.
   - (b) Tune GMRES tolerance/iterations so it converges (e.g. `gmres(Amult, RHS, 1000, tol, 1)`).
   - (c) Modify the right-hand side to account for the Dirichlet data on the outer box.
5. **Verify convergence:** the solution should still converge in Δx, but at a reduced
   order of **Δx^1.5** when the IB is included.
6. The result should look identical to the original solution (see Fig. 2 of the assignment).

### Problem 3 — An Almost Minimal Surface
Solve the same problem with f(x, y) = 0 and homogeneous Dirichlet boundary conditions
on the box, keeping the internal condition on the immersed boundary:

- ∇²u = S^ib(q), u = 0 on ∂ box, u(x_ib, y_ib) = cos(2x_ib) sin(2y_ib)

The solution's shape is closely related to a *minimal surface* — the shape a soap film
would take on the given wire-frame.

### Problem 4 — Charlie's Rule (play with the code!)
Pick something fun to show off, e.g.:

1. Speed up S^ib / J^ib by evaluating the Gaussian kernel only where it exceeds a tolerance.
2. Add multiple immersed shapes to the simulation.
3. Add pseudo-dynamics by moving/deforming the immersed boundary in time and re-solving
   (make videos).
4. Improve GMRES convergence with a good initial guess and/or a **preconditioner**
   (hint: Green's functions).
5. Extend the code to impose a **Neumann condition** ∂u/∂n on the immersed boundary
   instead of Dirichlet.

*(Bonus idea from the assignment footnote: the true minimal-surface PDE
`Div( ∇u / sqrt(1 + |∇u|²) ) = 0` is nonlinear but doable — potential project material.)*

## Files

- `PDES_IB_HW2.pdf` — original assignment handout (gitignored).
