It is 2D-ising model with square lattice.\
There are some advises for Classical Monte Carlo.
- To set a number of Markov chain to reduce the error. (mcset)
- Each Markov chain contains a list of measurements, and takes measurement for evergy Nx*Ny sites. (reduction correlation of measurements)
- For autocorrection time, we will take all observations of a Markov chain. However, only first finite interval is valid and those datas afterward are noisy.
- To select a nice updating scheme.

For updating scheme:
- standard simulation: local update
- Swendsen-Wang algorithm: cluster update

For Swenden-Wang algorithm(wiki):\
The problem of the critical slowing-down affecting local processes is of fundamental importance in the study of second-order phase transitions (like ferromagnetic transition in the Ising model), as increasing the size of the system in order to reduce finite-size effects has the disadvantage of requiring a far larger number of moves to reach thermal equilibrium. Indeed the correlation time \tau  usually increases as L^z with or greater; since, to be accurate, the simulation time must be larger than \tau, this is a major limitation in the size of the systems that can be studied through local algorithms. SW algorithm was the first to produce unusually small values for the dynamical critical exponents: 
- z=0.35 for the 2D Ising model (  z=2.125 for standard simulations)
- z=0.75 for the 3D Ising model, as opposed to z=2.0 for standard simulations.
