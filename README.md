# Diffusion Limited Aggregation (DLA) – Optimized C Implementation

High-performance simulation of Diffusion Limited Aggregation (DLA) written in C.  
The code implements several algorithmic optimizations that allow large off-lattice clusters to be simulated efficiently while collecting detailed statistical information about the growing interface.

---

# Overview

Diffusion Limited Aggregation is a stochastic growth process where particles perform random walks and stick when they touch an existing cluster.  
The resulting structures are fractal and appear in many physical systems such as electrodeposition, dielectric breakdown, mineral deposition, and viscous fingering.

This implementation focuses on performance and large-scale simulations.

Main goals of the code:

• Efficient random walker simulation  
• Fast collision detection  
• Reduced computational cost far from the cluster  
• Statistical analysis during growth  

---

# Main Optimizations

## Distance Grid (Omega)

A distance field stores an estimate of the distance between any position and the cluster.

This allows the walker to:

• take large steps when far from the cluster  
• reduce step size near the interface  
• avoid unnecessary collision checks

This significantly reduces runtime compared to naive DLA implementations.

---

# Spatial Hash Grid

Particles are stored in a coarse spatial grid.

Data structures used:

W_grid  
Maps spatial cells to particle blocks.

Nk_counts  
Number of particles stored in each cell.

F_indices  
Contiguous array containing particle indices.

Instead of checking collisions against all particles, the algorithm only checks nearby cells.

This reduces complexity from roughly:

O(N)

to approximately

O(local density)

---

# Continuous Collision Detection

When a walker is close to the cluster the code analytically computes the exact collision point.

The algorithm solves the quadratic equation describing the intersection between the walker trajectory and the particle.

Functions used:

solve_quadratic  
find_collision_distance

This improves both accuracy and performance.

---

# Grid System

Three grids are used in the simulation.

Y grid  
Stores particle positions on a lattice representation of the cluster.

Omega grid  
Stores the distance to the nearest particle.

Vicinity grid  
Precomputed distance template used to update the Omega grid efficiently.

This avoids recomputing distances globally.

---

# Interface Analysis

During the simulation the program measures properties of the interface:

• Radius distribution  
• Active growth zone  
• Maximum radius per angular sector  
• Roughness  
• Radius of gyration  

Statistics calculated include:

• Mean  
• Variance  
• Skewness  
• Kurtosis  
• Higher order moments

These are computed layer by layer as the cluster grows.

---

# Simulation Parameters

Main parameters defined in the code:

L = 8192  
Simulation box size

T = 100000  
Number of particles

layers = 100  
Number of measurement layers

D_MAX = 90  
Maximum reliable distance in the distance grid

Particle diameter:

2.0

Launch radius:

max_r + rmax

---

# Random Number Generator

The simulation uses a fast xorshift random generator.

Seed is generated from:

• system time  
• process ID  
• optional job ID

This allows running multiple simulations in parallel without correlations.

---

# Compilation

Compile using gcc:

gcc -mcmodel=large -Wall -O3 dla_opt.c -lm -o dla

---

# Running the Simulation

Run normally:

./dla

Run with job identifier (useful for clusters):

./dla job_id

---

# Output

The program can store:

Particle coordinates

x y

Radius distributions

r

These files can be used later for visualization or statistical analysis.

---

# Applications

• Fractal growth  
• Statistical physics  
• Monte Carlo simulations  
• Pattern formation  
• Complex systems research

---

# Author

Derlis Alfonso

Scientific computing  
C programming  
Monte Carlo simulations
