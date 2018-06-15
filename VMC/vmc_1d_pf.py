import numpy as np
from numpy import linalg
import copy, time, sys, os
import multiprocessing

epsilon = 1e-15

def generating_w(nx, params, bc):
    '''
    Constructe the mean field Hamiltonian and give w
    :param nx: size of x
    :param params: [t, ds, dt, d]
    :param bc: boundary condition = 0 (open) or 1 (periodic)
    :return: w
    '''
    N = nx
    t = 0
    ds, dt, d = params
    print('t , ds, dt, d = ', t, params)
    hkin = np.zeros((N,N)) # N = nx
    hsinglet = np.zeros((N,N))
    htriplet = np.zeros((N,N))
    ht = np.zeros((N,N))
    # Hamiltonian
    for px in range(nx):
        # position (px, px+1) by mod(): a % b
        a0 = px
        a1 = (px + 1) % (nx)

        hkin[a0, a1] = -t
        hkin[a1, a0] = -t

        hsinglet[a0, a1] = ds
        hsinglet[a1, a0] = ds

        htriplet[a0, a1] = -dt
        htriplet[a1, a0] = dt

        ht[a0, a1] = -d
        ht[a1, a0] = d

    h = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]), hkin)
    h += np.kron(np.array([[0, 0, 0, 1], [0, 0, -1, 0], [0, -1, 0, 0], [1, 0, 0, 0]]), hsinglet)
    h += np.kron(np.array([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]), htriplet)
    h += np.kron(np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]), ht)

    # h += np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), (1e-6) * np.random.random_sample((N, N)))

    ew, ev = np.linalg.eigh(h)
    u = ev[:(2 * N), :(2 * N)]
    v = -ev[(2 * N):, :(2 * N)]  # take u and v from same column
    u, v = np.conj(u.T), np.conj(v.T)
    ##
    P , sigma, QT = np.linalg.svd(u, full_matrices=True)  # becuase U is singular !!!!!
    ##
    print('check 2*N = ', np.trace(u.dot(np.conj(u.T)) + v.dot(np.conj(v.T))))
    print('check 0 = ', np.trace(u.dot(v.T) + v.dot(u.T)))
    ww = np.dot(np.dot(np.conj(QT.T), np.dot(np.linalg.pinv(np.diag(sigma)), np.conj(P.T))) , v) # np.dot(np.linalg.inv(u), v)
    print('sum W = ', ww.shape)
    check_w = np.max(np.abs(ww + ww.T))
    print('Max( Abs(W + W.T))) = ', check_w)
    print("Max (W) = ", np.max(ww))
    if( check_w >= 1e-8):
        print("W + W.T is not zero matrix!")
        check_w = False
    return ww, check_w

def vmc_mc(nx, ww):
    '''
    The body of MC
    :param nx:
    :param ww:
    :return:
    '''
    N, set_mc = nx, 10
    # random initial state = [position: 0, 1, 2, ...]
    state = np.arange(N)
    np.random.shuffle(state)
    print('State = ', state)
    print('Check random inital state (=(N-1)*N/2) = ', sum(state))
    print('Parity of state = ', perm_parity(state))
    proba = probability_pf(N, state, ww)
    print('Proba = ', proba)
    
    cores = multiprocessing.cpu_count()
    print('Total cores = ', cores)
    pool = multiprocessing.Pool(processes=cores)

    # for p1 in range(set_mc): # set_MC
    #    zz_result[p1, :, :], xy_result[p1, :, :] = parallel_mc(nx, state, ww, mc_step=1000)
    
    tasks = [(nx, state, ww, 500) for x in range(set_mc)]
    result = pool.starmap(parallel_mc, tasks)
    result = np.array(result) # list to array
    print('Size of Result of starmap = ', result.shape)
    pool.close()
    pool.join()
    
    zz_result, xy_result = result[:, 0, :, :], result[:, 1, :, :]

    print('ZZ size , XY size = ', zz_result.shape, xy_result.shape)
        
    mean_zz = zz_result.mean(axis=0)
    mean_xy = xy_result.mean(axis=0)
    std_zz = zz_result.std(axis=0)
    std_xy = xy_result.std(axis=0)
    print('Energy per spin <S_0 S_i >= ')
    print('Mean ZZ and Mean XY')
    print(mean_zz, '\n', mean_xy)
    print('Std ZZ and Std XY')
    print(std_zz, '\n', std_xy)
    mean = np.vstack((mean_zz, mean_xy))
    mean = np.vstack((mean, std_zz))
    mean = np.vstack((mean, std_xy))
    np.savetxt('zz_xy_'+str(nx)+'.csv', mean, delimiter=',')

    # position (px, px + 1) by mod(): a % b
    E_zz, E_xy = [], []
    for px in range(nx):
            a0 = px
            a1 = (px + 1) % (nx)
            E_zz.append(mean_zz[a0, a1])
            E_xy.append(mean_xy[a0, a1])
    print("N, mean E_zz, mean E_xy = ", N, np.mean(E_zz), np.mean(E_xy))
    print("N, std E_zz, std E_xy = ", N, np.std(E_zz), np.std(E_xy))
    return np.mean(E_zz), np.mean(E_xy), np.std(E_zz), np.std(E_xy)


def probability_pf(N, state, ww):
    alpha = np.zeros((N, N))
    for p1 in range(N):
        for p2 in range(N):
            if( p1 < (N/2)):
                if( p2 < (N/2)):
                    alpha[state[p1], state[p2]] = ww[state[p1], state[p2]] # up up
                else:
                    alpha[state[p1], state[p2]] = ww[state[p1], N + state[p2]] # up down
            else:
                if (p2 < (N / 2)):
                    alpha[state[p1], state[p2]] = ww[N + state[p1], state[p2]]  # down up
                else:
                    alpha[state[p1], state[p2]] = ww[N + state[p1], N + state[p2]]  # down down

    # print('Sum of alpha = ', sum(sum(alpha + alpha.T)))
    # print(abs((alpha + alpha.T).max()))

    # for singlet pairing
    # alpha = np.zeros((int(N/2), int(N/2)))
    # for p1 in range(int(N/2)):
    #     for p2 in range(int(N/2)):
    #         alpha[p1, p2] = ww[state[p1], N + state[p2 + int(N/2)]]

    # return perm_parity(state)*np.linalg.det(alpha)
    return pfaffian_LTL(alpha)

def perm_parity(lst):
    '''\
    Given a permutation of the digits 0..N in order as a list,
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i, len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity

def parallel_mc(nx, state, ww, mc_step):
    '''
    Main MC procedure and for parallel
    :param nx:
    :param state:
    :param ww:
    :param mc_step:
    :return:
    '''
    np.random.seed() # it is very important for parallel, seed are same if no this line
    N = nx
    xy, zz = np.zeros((N, N)), np.zeros((N, N))
    np.random.shuffle(state)
    # S_0 to S_p1
    for p0 in range(N):
        for p1 in range(N):
            txy, tzz = 0, 0
            proba = probability_pf(N, state, ww)
            # MC step
            for p2 in range(mc_step*N): # every N step, do measurement
                rup = np.random.randint(0, N/2) # Return random integers from low (inclusive) to high (exclusive), interval [low, high).
                rdown = np.random.randint(N/2, N)
                new_state = copy.copy(state)
                new_state[rup] = state[rdown]
                new_state[rdown] = state[rup]
                new_proba = probability_pf(N, new_state, ww)

                rand_t = np.random.random()
                # print(state, proba, new_proba, (abs(new_proba)/abs(proba))**2, rand_t)
                if (abs(new_proba) >= epsilon):
                    if( (abs(new_proba)/abs(proba))**2 >= rand_t):
                        state = new_state
                        proba = new_proba
                if( (p2%N) == 0):
                    test = measurement_zz(N, state, p0, p1)
                    tzz += test/4.0
                    txy += 0.5*(measurement_xy(N, state, ww, p0, p1, test)/proba)
            # done for mc_step
            zz[p0, p1] += tzz/mc_step
            xy[p0, p1] += txy/mc_step
    return zz, xy


def measurement_zz(N, state, position1, position2):
    '''
    Calculate < S_0 S_{position2} >
    :param nx:
    :param ny:
    :param state:
    :param position2:
    :return:
    '''
    test = 1
    for p1 in range(int(N/2)):
        if( state[p1] == position1):
            test *= -1
            break
    for p1 in range(int(N / 2)):
        if (state[p1] == position2):
            test *= -1
            break
    return test

def measurement_xy(N, state, ww, position1, position2, test_from_zz):
    if( test_from_zz > 0):
        return 0
    xy_state = copy.copy(state)
    for p1 in range(N):
        if (state[p1] == position1):
            aa = p1
        if (state[p1] == position2):
            bb = p1
    temp = xy_state[aa]
    xy_state[aa] = xy_state[bb]
    xy_state[bb] = temp

    xy_proba = probability_pf(N, xy_state, ww)
    if (abs(xy_proba) <= epsilon):
        xy_proba = epsilon
    return xy_proba

def pfaffian_LTL(A, overwrite_a=False):
    """ pfaffian_LTL(A, overwrite_a=False)

    Compute the Pfaffian of a real or complex skew-symmetric
    matrix A (A=-A^T). If overwrite_a=True, the matrix A
    is overwritten in the process. This function uses
    the Parlett-Reid algorithm.
    """
    #Check if matrix is square
    assert A.shape[0] == A.shape[1] > 0
    #Check if it's skew-symmetric
    assert np.max(np.abs(A+A.T)) < 1e-10 #  abs((A+A.T).max()) < 1e-10

    n = A.shape[0]
    A = np.asarray(A)  #the slice views work only properly for arrays

    #Quick return if possible
    if n%2==1:
        return 0

    if not overwrite_a:
        A = A.copy()

    pfaffian_val = 1.0

    for k in range(0, n-1, 2):
        #First, find the largest entry in A[k+1:,k] and
        #permute it to A[k+1,k]
        kp = k+1+np.abs(A[k+1:,k]).argmax()

        #Check if we need to pivot
        if kp != k+1:
            #interchange rows k+1 and kp
            temp = A[k+1,k:].copy()
            A[k+1,k:] = A[kp,k:]
            A[kp,k:] = temp

            #Then interchange columns k+1 and kp
            temp = A[k:,k+1].copy()
            A[k:,k+1] = A[k:,kp]
            A[k:,kp] = temp

            #every interchange corresponds to a "-" in det(P)
            pfaffian_val *= -1

        #Now form the Gauss vector
        if A[k+1,k] != 0.0:
            tau = A[k,k+2:].copy()
            tau /= A[k,k+1]

            pfaffian_val *= A[k,k+1]

            if k+2<n:
                #Update the matrix block A(k+2:,k+2)
                A[k+2:,k+2:] += np.outer(tau, A[k+2:,k+1])
                A[k+2:,k+2:] -= np.outer(A[k+2:,k+1], tau)
        else:
            #if we encounter a zero on the super/subdiagonal, the
            #Pfaffian is 0
            return 0.0

    return pfaffian_val


if __name__ == '__main__':
    print('=' * 50)
    print('Start with ', time.asctime(time.localtime(time.time())))
    t1 = time.time()
    print('=' * 50)

    nx = 6
    N = nx
    data_result_zz, data_result_xy, data_result_std_zz, data_result_std_xy = [], [], [], []
    print('The system size N = ', N, '(', nx, ' )')
    x = np.linspace(0, 2, 21)
    for p1 in [0.001]:
        print('=' * 50)
        ww, check_w = generating_w(nx, [p1, 1e-6, 1e-7], 1) # params: t=-1 [ds, dt, d]
        if(check_w == False):
            data_result_zz.append(0.0)
            data_result_xy.append(0.0)
            data_result_std_zz.append(0.0)
            data_result_std_xy.append(0.0)
            continue 
        zz, xy, stdzz, stdxy = vmc_mc(nx, ww)
        data_result_zz.append(zz)
        data_result_xy.append(xy)
        data_result_std_zz.append(stdzz)
        data_result_std_xy.append(stdxy)
    print("Results: ", data_result_zz)
    print("Results: ", data_result_xy)
    print("Results: ", data_result_std_zz)
    print("Results: ", data_result_std_xy)
    
    #--- save data ----
    f = open('nx_%s_paras_data.csv'%(nx), 'wb')
    np.savetxt(f, x[np.newaxis], delimiter=',') # if 1-D, a[np.newaxis]
    np.savetxt(f, np.asarray(data_result_zz)[np.newaxis], delimiter=',') # if 1-D, a[np.newaxis]
    np.savetxt(f, np.asarray(data_result_xy)[np.newaxis], delimiter=',') # if 1-D, a[np.newaxis]
    np.savetxt(f, np.asarray(data_result_std_zz)[np.newaxis], delimiter=',') # if 1-D, a[np.newaxis]
    np.savetxt(f, np.asarray(data_result_std_xy)[np.newaxis], delimiter=',') # if 1-D, a[np.newaxis]
    f.close()

    print('=' * 50)
    t2 = time.time()
    print('time = ', t2 - t1, ' s ')
    print('time = ', (t2 - t1) / 60, ' mins ')
    print('time = ', (t2 - t1) / (60 * 60), 'hours')
    print('time = ', (t2 - t1) / (60 * 60 * 24), ' days')
    print('End with ', time.asctime(time.localtime(time.time())))
    print('=' * 50)
