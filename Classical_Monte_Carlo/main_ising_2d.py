import numpy as np
import time, os, copy
import multiprocessing
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["savefig.directory"]=os.chdir(os.getcwd())
mpl.rcParams.update({'font.size': 22})
import ising_2d_v2

def testing():
    # run_mc(self, T, mcset, total_mc_step, update = 'standard', Do_Autocorrelation_time = False)
    # (nx, ny) = (8, 8) AND -1 for ferro, +1 for anti-ferro
    # check form T = 1.5 to 4
    n = 4
    nx, ny, J = n, n, -1
    test_model = ising_2d_v2.ising_model(nx, ny, J)
    print(test_model.__doc__)
    m = 30
    updating = 'standard'
    t_list = np.linspace(1.5, 4, m)
    
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    print("CPU cores = ", cores)
    tasks = [(t_list[x], cores, 1000, updating) for x in range(m)]
    data_out = pool.starmap(test_model.run_mc, tasks)
    pool.close()
    pool.join()
    data_out = np.asarray(data_out)
    print("Data output size = ", data_out.shape)
    
    np.savetxt('data_ising_2d_{0}_{1}_{2}_{3}.csv'.format(nx, ny, J, updating), data_out, delimiter=',')
    

def testing_auto_corr_time():
    # run_mc(self, T, mcset, total_mc_step, update = 'standard', Do_Autocorrelation_time = False)
    # (nx, ny) AND -1 for ferro, +1 for anti-ferro
    # check form T = 1.5 to 4
    n = 8
    nx, ny, J = n, n, -1
    test_model = ising_2d_v2.ising_model(nx, ny, J)
    print(test_model.__doc__)
    m = 1
    total_mc = 1000
    updating = 'standard'

    data1, data2 = np.zeros((m, total_mc)), np.zeros((m, total_mc))
    for p1 in range(m):
        data1[p1, :], data2[p1, :]= test_model.run_mc(T = 2.65, mcset = 1, total_mc_step = total_mc, update=updating, Do_Autocorrelation_time = True)
    
    data1 = np.mean(data1, axis=0)
    data2 = np.mean(data2, axis=0)

    np.savetxt('data_2d_auto_corr_time_{0}_{1}_{2}_{3}.csv'.format(nx, ny, J, updating), np.vstack((data1, data2)), delimiter=',')
            
    
if __name__ == '__main__':
    print(time.asctime(time.localtime(time.time())))
    print("="*50)
    t1 = time.time()
    
    # testing()
    testing_auto_corr_time()
    
    print("="*50)
    t2 = time.time()
    print('time = ', t2-t1, ' s ')
    print('time = ', (t2-t1)/60, ' mins ')
    print(time.asctime(time.localtime(time.time())))
    print("="*50)
