It is 2D-ising model with square lattice.\
There are some advises for Classical Monte Carlo.
- To set a number of Markov chain to reduce the error. (mcset) (Q: why and how?) https://web.northeastern.edu/afeiguin/phys5870/phys5870/node71.html
- Each Markov chain contains a list of measurements, and takes measurement for evergy Nx*Ny sites. (reduction correlation of measurements)
- For autocorrection time, we will take all observations of a Markov chain. However, only first finite interval is valid and those datas afterward are noisy. (Q: why?)
- To select a nice updating scheme.
- When you want to study the model clearly, you should save all measurements of all Markov Chains. It is quite time-consuming for large size, so just save all information and then do analysis later.

For updating scheme:
- standard simulation: local update
- Swendsen-Wang algorithm: cluster update

For Swenden-Wang algorithm(wiki) and the advance version -- Wolff algorithm :\
The problem of the critical slowing-down affecting local processes is of fundamental importance in the study of second-order phase transitions (like ferromagnetic transition in the Ising model), as increasing the size of the system in order to reduce finite-size effects has the disadvantage of requiring a far larger number of moves to reach thermal equilibrium. Indeed the correlation time \tau  usually increases as L^z with or greater; since, to be accurate, the simulation time must be larger than \tau, this is a major limitation in the size of the systems that can be studied through local algorithms. SW algorithm was the first to produce unusually small values for the dynamical critical exponents: 
- z=0.35 for the 2D Ising model (  z=2.125 for standard simulations)
- z=0.75 for the 3D Ising model, as opposed to z=2.0 for standard simulations.

For Wolff algorithm:\
Let's the spin Hamiltonian as
- H = -J \sum_{i,j} S_i S_j, J>0 and ferro case
- To select randomly a site, add it into a cluster C.
- To find all neighbors of this site except on cluster C, activate it by same spin sign and P = max{0, 1-exp(-2 \beta J}, for ferror case.
- If activate,put it into cluster C.
- To find out all possible neighbors
- Flip all spins in the cluster
- Tips: set of cluster C includes all possible sites and set of Remaining includes C but not calculated yet.

In setting P = max{0, 1-exp(-2 \beta}, the accpectance rate becomes 1.(However, Wolff will take more computation time...)
- under same conditions, standard takes 1.2 mins and wolff takes 58 mins
- of course, the autocorrelation time will show that wolff can take less steps to give same accuracy.
