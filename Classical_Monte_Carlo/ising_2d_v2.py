import numpy as np
import time, os, copy
import multiprocessing
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["savefig.directory"]=os.chdir(os.getcwd())
mpl.rcParams.update({'font.size': 22})

class ising_model:
    """
    The 2D ising model
    (nx, ny) = dimension
    J = -1 (ferro) or 1 (anti-ferro)
    
    model.run_mc(self, T, mcset, total_mc_step, update = 'standard' or 'wolff', Do_Autocorrelation_time = True or False)
    
    Energy <E>, Magnetization <|M|>, <E^2>, <|M|^2>
    
    output = [mean_energy, std_energy, mean_m, std_m, cv, std_cv, chi_m, std_chim]
    """
    def __init__(self, nx, ny, J):
        self.nx = nx
        self.ny = ny
        self.J = J 
    
    def run_mc(self, T, mcset, total_mc_step, update = 'standard', Do_Autocorrelation_time = False):
        beta = 1.0/T
        finial_e, finial_m = np.zeros(mcset), np.zeros(mcset)
        finial_e2, finial_m2 = np.zeros(mcset), np.zeros(mcset)
        """
        There are [mcset] Markov chains, and we takes those average as \sum_i=mcset <E>_i
        It is important for reducing the std or errors
        """
        print("="*50)
        print(" Starting MC ... ")
        print("(Nx, Ny) = ({0}, {1}), T = {2}, mcset = {3}, total_mc_step = {4}, update = {5}, do_Autocorrelation_time = {6}, J = {7}".\
              format(self.nx, self.ny, T, mcset, total_mc_step, update, Do_Autocorrelation_time, self.J))
        
        for p2 in range(mcset):
            # initialize spin configuration
            s = np.random.choice([-1, 1], size=(self.nx, self.ny))
            energy, magnetization = [], []
            
            # Warm up
            for __ in range(100*self.nx*self.ny):
                s = self.updating_scheme(s, beta, update) # 'wolff'

            # Markov-chain starts
            # for standard update, take meansurement for evergy nx*ny
            if update[0] == 's':
                for p1 in range(total_mc_step*self.nx*self.ny):
                    s = self.updating_scheme(s, beta, update)
                    # Measurement
                    if (np.mod(p1, self.nx*self.ny) == 0):
                        energy.append(self.cal_energy(s)/(self.nx*self.ny))
                        magnetization.append(self.cal_m(s)/(self.nx*self.ny))      
            # for wolff update, take meansurement for evergy mc single step, (btw up to you la)
            elif update[0] == 'w':
                for p1 in range(total_mc_step):
                    s = self.updating_scheme(s, beta, update)
                    # Measurement
                    energy.append(self.cal_energy(s)/(self.nx*self.ny))
                    magnetization.append(self.cal_m(s)/(self.nx*self.ny))
            
            finial_e[p2] = np.mean(energy)
            finial_m[p2] = np.mean(magnetization)
            finial_e2[p2] = np.mean(np.asarray(energy)**2)
            finial_m2[p2] = np.mean(np.asarray(magnetization)**2)
        
        # Autocorrelation time, take the last Markov-chain, only one Markov-chain!
        # True or False
        if( Do_Autocorrelation_time == True): # True or False
            ac_energy = self.auto_corr_time(energy, np.asarray(energy)**2)
            ac_m = self.auto_corr_time(magnetization, np.asarray(magnetization)**2)
            return ac_energy, ac_m
        else:
            mean_energy, std_energy = np.mean(finial_e), np.std(finial_e)
            mean_m, std_m = np.mean(finial_m), np.std(finial_m)
            cv, std_cv = self.cal_variance(finial_e, finial_e2)/T**2
            chi_m, std_chim = self.cal_variance(finial_m, finial_m2)/T
            return np.array([mean_energy, std_energy, mean_m, std_m, cv, std_cv, chi_m, std_chim])
    
    def updating_scheme(self, s, beta, updating_method):
        # standard method
        if updating_method[0] == 's': 
            x1, x2 = np.random.randint(self.nx), np.random.randint(self.ny)
            s_change = s[x1, x2]
            # let's set it as square lattice with PBC
            neighbr = s[(x1 + 1) % self.nx, x2] + s[x1, (x2 + 1) % self.ny] + s[(x1 - 1) % self.nx, x2] + s[x1, (x2 - 1) % self.ny]
            # new_E - old_E
            cost_energy = -2*self.J*s_change*neighbr
            if rand() < np.exp(-cost_energy*beta):
                s[x1, x2] *= -1
            return s
        # Wolff algorithm
        elif updating_method[0] == 'w':
            new_s = copy.copy(s)
            cluster_list, remaining = [], []
            x1, x2 = np.random.randint(self.nx), np.random.randint(self.ny)
            new_s[x1, x2] *= -1
            cluster_list.append([x1, x2])
            remaining.append([x1, x2])
            nxy = np.zeros((4,2))
            while(len(remaining) != 0):
                # to take the first position of remaining
                r1, r2 = remaining[0][0], remaining[0][1]
                nxy[0, :] = np.array([(r1 + 1) % self.nx, r2])
                nxy[1, :] = np.array([r1, (r2 + 1) % self.ny])
                nxy[2, :] = np.array([(r1 - 1) % self.nx, r2])
                nxy[3, :] = np.array([r1, (r2 - 1) % self.ny])
                # sign of spin must be same with s[x1, x2]
                for p1 in range(4):
                    if(s[x1, x2] == s[int(nxy[p1, 0]), int(nxy[p1, 1])]) and (not([nxy[p1, 0], nxy[p1, 1]] in cluster_list)):
                        if rand() < (1-np.exp(2*self.J*beta)): # if ferro, it is 1 - exp(-2*beta)
                            new_s[int(nxy[p1, 0]), int(nxy[p1, 1])] *= -1
                            cluster_list.append([nxy[p1, 0], nxy[p1, 1]])
                            remaining.append([nxy[p1, 0], nxy[p1, 1]])
                remaining.remove([r1, r2])    
            return new_s

    
    def cal_energy(self, s):
        temp = 0.0
        for x1 in range(self.nx):
            for x2 in range(self.ny):
                temp += (s[(x1 + 1) % self.nx, x2] + s[x1, (x2 + 1) % self.ny])*s[x1, x2]
        return self.J*temp
    
    def cal_m(self, s):
        # <|M|> measurement
        if (self.J > 0):
            # if it is Anti-Ferromagnet
            temp = 0.0
            for x1 in range(self.nx):
                for x2 in range(self.ny):
                    temp += (-1)**(x1+x2) * s[x1, x2]
            return np.abs(temp)
        else:
            # if it is Ferromagnet
            return np.abs(np.sum(s))
        
    
    def cal_variance(self, data_list, data_list_2):
        # like heat-capacity and spin susceptibility
        # be carefull the equration !! and over the T**2 or T
        output = data_list_2 - data_list**2
        return np.mean(output), np.std(output)
    
    def auto_corr_time(self, data_list, data_list_2):
        # autocorrelation time from 0 (self) to N
        data_list = np.asarray(data_list)
        auto_time, mean_data = [], np.mean(data_list)
        length_data = data_list.shape[0]
        # p1 from 0 to 10 or 20, it can't be large number. Otherwise, it makes wrong.
        for p1 in range(length_data):
            temp = 0
            for p2 in range(length_data - p1):
                temp += data_list[p2]*data_list[p2 + p1]
            temp = temp/(length_data - p1) - mean_data**2
            auto_time.append(temp)
        return np.asarray(auto_time)/(np.mean(data_list_2) - mean_data**2)
        