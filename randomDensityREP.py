from multiprocessing import Process, Array
import numpy as np
from ctypes import c_float
from sys import argv, maxsize
from subprocess import check_output

np.set_printoptions(threshold=maxsize)
import re
import threading

k = int(argv[1])
L = int(argv[2])
SIZE = 1000000


def main():
    if argv[4] == 'i':
        outputs = densityI(seqpath=argv[5], reps=int(argv[3]))
    else:
        REP = int(argv[3])
        outputs = Array(c_float, np.empty([REP]))
        # densitystart = time()
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
    exit(0)


def densityR(OPtable, idx):
    np.random.seed(idx)
    s = np.array2string(np.random.randint(low=0, high=4, size=SIZE, dtype='i1'),
                        separator="", suffix="", prefix="")[1:-1].replace('\n', '').replace(' ', '')
    currmzeridx = findNewOne(s, 0)
    currmzer = s[currmzeridx: currmzeridx + k]
    n = 1
    lidx = 1
    while lidx < SIZE - L + 1:
        if lidx - 1 == currmzeridx:
            currmzeridx = findNewOne(s, lidx)
            currmzer = s[currmzeridx: currmzeridx + k]
            n += 1
        else:
            newkmeridx = lidx + L - k
            newkmer = s[newkmeridx: newkmeridx + k]
            currvalue = hash(currmzer)
            newvalue = hash(newkmer)
            if currvalue > newvalue:
                currmzeridx = newkmeridx
                currmzer = newkmer
                n += 1
        lidx += 1
    dens = n / (SIZE - k + 1)
    OPtable[idx] = dens


def densityI(seqpath, reps):
    outputs = np.empty([reps])
    for i in range(reps):
        cmd = 'python3 randomF.py ' + str(k) + ' ' + str(L) + ' ' + seqpath
        output = check_output(cmd, shell=True)
        outputs[i] = float(output)
    return outputs


# def density(size, s, OPtable, idx):
#     if s is None:
#         s = np.array2string(np.random.randint(low=0, high=4, size=size, dtype='i1'), separator="", suffix="",
#                             prefix="")[1:-1].replace('\n', '').replace(' ', '')
#     currmzeridx = findNewOne(s, 0)
#     currmzer = s[currmzeridx: currmzeridx + k]
#     lidx = 0
#     n = 1
#     lidx += 1
#     while lidx < size - L + 1:
#         # print("lidx is ", lidx)
#         if lidx - 1 == currmzeridx:
#             currmzeridx = findNewOne(s, lidx)
#             currmzer = s[currmzeridx: currmzeridx + k]
#             n += 1
#         else:
#             newkmeridx = lidx + L - k
#             newkmer = s[newkmeridx: newkmeridx + k]
#             currvalue = hash(currmzer)
#             newvalue = hash(newkmer)
#             if currvalue > newvalue:
#                 currmzeridx = newkmeridx
#                 currmzer = newkmer
#                 n += 1
#         lidx += 1
#     OPtable[idx] = n / (size - k + 1)


def findNewOne(s, idx):
    lmer = s[idx: idx + L]
    bestKmerSoFaridx = 0
    bestValueSoFar = hash(lmer[:k])
    i = 1
    while i < L - k + 1:
        currKmer = lmer[i: i + k]
        currValue = hash(currKmer)
        if currValue < bestValueSoFar:
            bestValueSoFar = currValue
            bestKmerSoFaridx = i
        i += 1
    return idx + bestKmerSoFaridx


if __name__ == '__main__':
    main()
