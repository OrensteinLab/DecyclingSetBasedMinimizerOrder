from multiprocessing import Process, Array
import numpy as np
from ctypes import c_float
from sys import argv, maxsize
from subprocess import check_output
import functions
from math import sin
# from time import time

np.set_printoptions(threshold=maxsize)

k = int(argv[1])
L = int(argv[2])
SIZE = 1000000


def main():
    if argv[4] == 'i':
        outputs = densityI(seqpath=argv[5], reps=int(argv[3]))
    else:
        global coeffs
        u = 3.141592653589793 * 2.0 / k
        coeffs = [sin(i * u) for i in range(1, k)]
        REP = int(argv[3])
        outputs = Array(c_float, np.empty([REP]))
        if REP == 1:
            densityR(outputs, 0)
            print(f'{outputs[0]} 0')
            exit(0)
        else:
            processes = []
            for i in range(REP):
                p = Process(target=densityR, args=(outputs, i))
                p.start()
                processes.append(p)
            for process in processes:
                process.join()
    print(f'{np.mean(outputs)} {np.std(outputs)}')
    return 0


def densityR(OPtable, idx):
    np.random.seed(idx)
    s = np.array2string(np.random.randint(low=0, high=4, size=SIZE, dtype='i1'),
                        separator="", suffix="", prefix="")[1:-1].replace('\n', '').replace(' ', '')
    currmzeridx = functions.decfindNewOne(s, 0, k, L, compare, coeffs)
    currmzer = s[currmzeridx: currmzeridx + k]
    n = 1
    lidx = 1
    while lidx < SIZE - L + 1:
        if lidx - 1 == currmzeridx:
            currmzeridx = functions.decfindNewOne(s, lidx, k, L, compare, coeffs)
            currmzer = s[currmzeridx: currmzeridx + k]
            n += 1
        else:
            newkmeridx = lidx + L - k
            newkmer = s[newkmeridx: newkmeridx + k]
            if compare(newkmer, currmzer, coeffs) > 0:
                currmzeridx = newkmeridx
                currmzer = newkmer
                n += 1
        lidx += 1
    dens = n / (SIZE - k + 1)
    OPtable[idx] = dens



def densityI(seqpath, reps):
    outputs = np.empty([reps])
    for i in range(reps):
        cmd = 'python3 decmodF.py ' + str(k) + ' ' + str(L) + ' ' + seqpath
        output = check_output(cmd, shell=True)
        outputs[i] = float(output)
    return outputs



# def densityI(size, s, OPtable, idx):
#     print(f'idx is {idx}')
#     currmzeridx = functions.decfindNewOne(s, 0, k, L, compare, coeffs)
#     currmzer = s[currmzeridx: currmzeridx + k]
#     n = 1
#     lidx = 1
#     while lidx < size - L + 1:
#         if lidx - 1 == currmzeridx:
#             currmzeridx = functions.decfindNewOne(s, lidx, k, L, compare, coeffs)
#             currmzer = s[currmzeridx: currmzeridx + k]
#             n += 1
#         else:
#             newkmeridx = lidx + L - k
#             newkmer = s[newkmeridx: newkmeridx + k]
#             if compare(newkmer, currmzer) > 0:
#                 currmzeridx = newkmeridx
#                 currmzer = newkmer
#                 n += 1
#         lidx += 1
#     dens = n / (size - k + 1)
#     OPtable[idx] = dens


def compare(akmer, bkmer, coeffs):
    membership = [functions.checkIndecmod(akmer, coeffs), functions.checkIndecmod(bkmer, coeffs)]
    if membership[0] == membership[1]:
        avalue = hash(akmer)
        bvalue = hash(bkmer)
        if avalue == bvalue:
            return 0
        if avalue > bvalue:
            return 1
        else:
            return -1
    elif membership[0]:
        return 1
    return -1


if __name__ == '__main__':
    main()
