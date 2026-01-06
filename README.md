Project: 2D Poisson Solver for Discrete Charge Distribution
Overview

This project is a numerical simulation tool developed in MATLAB to solve the 2D Poisson equation for electrostatics:
∇²V = -ρ/ε

The solver computes the electric potential (V) and the electric field (E) on a rectangular domain based on user-defined point charges and boundary conditions. It features an interactive GUI for charge placement and real-time visualization of the field evolution.

Mathematical Foundation
1. The Poisson Equation

The project uses the Finite Difference Method (FDM) to solve the Poisson equation. On a uniform grid with spacing 'h', the continuous Laplacian operator is replaced by a 5-point stencil.

2. Numerical Methods (Iterative Solvers)

The user can choose between three iterative methods:

Jacobi Method: Simplest method, updates nodes using only values from the previous iteration. It is the slowest to converge.

Gauss-Seidel Method: Faster than Jacobi. It uses the most recently calculated values immediately within the same iteration.

Successive Over-Relaxation (SOR): The most efficient method. It uses a relaxation factor (ω) between 1 and 2 to accelerate convergence.

Recommendation: SOR is the preferred method for this project as it reduces the number of iterations required by nearly 90% compared to Jacobi.

Features
Interactive Charge Editor

Users can click on the domain to place point charges.

Each click prompts for a charge magnitude (q).

Charges are mapped to grid nodes using the relation: ρ = q / h².

Flexible Boundary Conditions

Dirichlet: Fixed potential (V) at the edge. Useful for modeling grounded walls or battery terminals.

Neumann: Fixed normal derivative (dV/dn) at the edge. Useful for modeling insulators or symmetry planes (Zero Flux).

Real-time Visualization

Potential Distribution: Heatmap/Contour plots of V(x,y).

Electric Field: Vector quiver plots showing E = -∇V.

Convergence Monitoring: Log-scale plot of the residual error vs. iteration count.

Algorithm Workflow

Algorithm 1: Grid Setup
Initialize the potential matrix (V) and charge density matrix (ρ). Map discrete point charges to the nearest grid nodes.

Algorithm 2: Iterative Solver
The core loop that updates Neumann boundaries, calculates interior node potentials using the selected method (Jacobi/GS/SOR), re-imposes Dirichlet conditions, and checks for convergence against a tolerance (τ).

Algorithm 3: Electric Field Computation
After the potential converges, the electric field components (Ex, Ey) are calculated using centered differences for interior nodes and one-sided differences for boundaries.

Algorithm 4: GUI Workflow
Handles user inputs, coordinate mapping for mouse clicks, and triggers the solver and visualization routines.

Installation and Usage

Open MATLAB.

Ensure all project files (.m files) are in the same folder.

Run the main script: PoissonSolver.m.

Follow the command window prompts or GUI fields to define Lx, Ly, Nx, Ny, and ε.

Click on the plot to place charges; right-click to finish.

Select the solver (SOR recommended) and define boundary conditions.

Click "Run" to observe the potential field converge in real-time.

Requirements

MATLAB (R2020a or later recommended).

Image Processing Toolbox (optional, for certain visualization enhancements).

Project Deliverables (Required)

Potential field heatmap/contour plot.

Electric field vector (quiver) plot.

Charge location overlay.

Convergence monitoring plot (Residual vs. Iteration).

Verification of boundary conditions.

Author: Nazmul Ibn Huda
Project Date: 2026
Course: Numerical Methods for Electromagnetics
