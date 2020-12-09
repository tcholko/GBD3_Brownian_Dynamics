# GBD3_Brownian_Dynamics

A sample of the source code from my PhD lab's (Chang, UCR) in-house Brownian dynamics programs, GeomBD3, on which I have done extensive work.
This repository contains only the source code files on which I've made my most significant contributions.

The features I've contributed to the program include:
- Spherical, rectangular or hexagonal simulation cell boundaries
- Periodic boundary conditions
- Planar, single-point or random-point initial ligand positioning
- Ability to use multiple ligand-receptor binding conditions which must be met either separately or simultaneously
- A feature to save all ligand coordinates and trajectory timers, allowing simulations to be stopped and restarted
  at any time. This makes it much more convenient to achieve converged simulation results in cases where many hundreds 
  of hours of wall-time on high-performace clusters may be needed
- Improved electrostatic and Lennard-Jones grid shape calculation
- Several trajectory analysis scripts e.g., diffusion coefficients, non-bonded interaction energy calculation, binding pathway
  analysis, ligand residence times, etc.
