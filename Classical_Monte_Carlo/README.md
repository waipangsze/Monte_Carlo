It is 2D-ising model with square lattice.\
There are some advises for Classical Monte Carlo.
- To set a number of Markov chain to reduce the error. (mcset)
- Each Markov chain contains a list of observations, and take observation for evergy Nx*Ny sites. (reduction correlation of observations)
- For autocorrection time, we will take all observations of a Markov chain. However, only first finite interval is valid and those datas afterward are noisy.
- To select a nice updating scheme. 
