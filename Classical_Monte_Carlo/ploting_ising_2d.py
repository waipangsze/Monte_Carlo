import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["savefig.directory"]=os.chdir(os.getcwd())
mpl.rcParams.update({'font.size': 22})

def ploting_mc(n, filename):
    
    data_out = np.loadtxt(open(filename, "r+b"), delimiter=",")
    plt.subplot(2,2,1)
    plt.plot(t_list, data_out[:, 0], 'bo--')
    plt.errorbar(t_list, data_out[:, 0], yerr=data_out[:, 1], fmt='o')
    plt.title("Energy with {0} x {1}".format(n, n))
    plt.grid()
    plt.subplot(2,2,2)
    plt.plot(t_list, data_out[:, 2], 'bo--')
    plt.errorbar(t_list, data_out[:, 2], yerr=data_out[:, 3], fmt='o')
    plt.title("<|M|>")
    plt.grid()
    plt.subplot(2,2,3)
    plt.plot(t_list, data_out[:, 4], 'bo--')
    plt.errorbar(t_list, data_out[:, 4], yerr=data_out[:, 5], fmt='o')
    plt.title("C_v")
    plt.grid()
    plt.subplot(2,2,4)
    plt.plot(t_list, data_out[:, 6], 'bo--')
    plt.errorbar(t_list, data_out[:, 6], yerr=data_out[:, 7], fmt='o')
    plt.title("chi_m")
    plt.grid()
    plt.show()
    plt.close()


def ploting_auto_corr_time(n, filename, cutoff = 20):

    data_out = np.loadtxt(open(filename, "r+b"), delimiter=",")
    sx, sy = data_out.shape
    data1, data2 = data_out[0, :], data_out[1, :]
    print("size = ", data_out.shape, data1.shape, data2.shape)
        
    """
    let's set a cutoff
    """
    plt.subplot(1,2,1)
    plt.plot(data1[:cutoff], 'bo--')
    plt.title("Energy with {0} x {1}".format(n, n))
    plt.grid()
    plt.subplot(1,2,2)
    plt.plot(data2[:cutoff], 'bo--')
    plt.title("M")
    plt.grid()
    plt.show()
    plt.close()
    
if __name__ == '__main__':
    n = 8
    filename = 'data_ising_2d_4_4_-1_standard.csv'
    # ploting_mc(n, filename)
    filename = 'data_2d_auto_corr_time_8_8_-1_standard.csv'
    ploting_auto_corr_time(n, filename, cutoff= 20)

    
    
