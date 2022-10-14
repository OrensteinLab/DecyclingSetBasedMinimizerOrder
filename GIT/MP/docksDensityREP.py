import numpy as np
from sys import argv, maxsize
from os.path import exists
from math import sin
from ctypes import c_float, c_uint32
from multiprocessing import Array, Process
import functions
from subprocess import check_output
from bitarray import bitarray
np.set_printoptions(threshold=maxsize)


k = int(argv[1])
L = int(argv[2])
SIZE = 1000000


def main():

    uhspath = "../docksdensitysets/PASHA" + str(k) + "_" + str(L) + ".txt"
    if not exists(uhspath):
        print('0 0')
        return 0
    if argv[4] == 'i':
        outputs = densityI(seqpath=argv[5], reps=int(argv[3]))
    else:
        REP = int(argv[3])
        global coeffs
        u = 3.141592653589793 * 2.0 / k
        coeffs = [sin(i * u) for i in range(1, k)]
        uhsset = bitarray('0' * pow(4, k))
        outputs = Array(c_float, np.empty([REP]))
        with open(uhspath, "r") as uhsfile:
            lines = uhsfile.readlines()
            for line in lines:
                transtable = line.maketrans("ACGT", "0123")
                translatedline = line.translate(transtable).replace('\n', '')
                idx = functions.toDeci(translatedline)
                uhsset[idx] = 1
        uhsarr = Array(c_uint32, functions.bitstoArray(k, uhsset))
        if REP == 1:
            densityR(outputs, 0, uhsarr)
            print(f'{outputs[0]} 0')
            exit(0)
        else:
            processes = []
            for i in range(REP):
                p = Process(target=densityR, args=(outputs, i, uhsset))
                p.start()
                processes.append(p)
            for process in processes:
                process.join()
    print(f'{np.mean(outputs)} {np.std(outputs)}')


def densityI(seqpath, reps):
    outputs = np.empty([reps])
    for i in range(reps):
        cmd = 'python3 docksF.py ' + str(k) + ' ' + str(L) + ' ' + seqpath
        output = check_output(cmd, shell=True)
        outputs[i] = float(output)
    return outputs



def densityR(OPtable, idx, uhsset):
    np.random.seed(idx)
    s = np.array2string(np.random.randint(low=0, high=4, size=SIZE, dtype='i1'),
                        separator="", suffix="", prefix="")[1:-1].replace('\n', '').replace(' ', '')
    currmzeridx = functions.uhsfindNewOne(s, 0, k, L, compare, coeffs, uhsset)
    currmzer = s[currmzeridx: currmzeridx + k]
    n = 1
    lidx = 1
    while lidx < SIZE - L + 1:
        if lidx - 1 == currmzeridx:
            currmzeridx = functions.uhsfindNewOne(s, lidx, k, L, compare, coeffs, uhsset)
            currmzer = s[currmzeridx: currmzeridx + k]
            n += 1
        else:
            newkmeridx = lidx + L - k
            newkmer = s[newkmeridx: newkmeridx + k]
            if compare(newkmer, currmzer, uhsset, coeffs) > 0:
                currmzeridx = newkmeridx
                currmzer = newkmer
                n += 1
        lidx += 1
    dens = n / (SIZE - k + 1)
    OPtable[idx] = dens


def compare(akmer, bkmer, uhsset, coeffs):
    membership = [functions.checkIndec(akmer, coeffs) or functions.checkmembership(akmer, uhsset), functions.checkIndec(bkmer, coeffs) or functions.checkmembership(bkmer, uhsset)]
    # print('akmer is ', akmer, 'and its in decycling:', membership[0], ' bkmer is ', bkmer, 'and its in decycling:', membership[1])
    if membership[0] == membership[1]:
        avalue = hash(akmer)
        bvalue = hash(bkmer)
        if avalue == bvalue:
            # print('theyre equal')
            return 0
        if avalue > bvalue:
            # print(akmer, ' is bigger')
            return 1
        else:
            # print(bkmer, ' is bigger')
            return -1
    elif membership[0]:
        # print('akmer', akmer, ' is in decset and bkmer', bkmer, ' is not')
        return 1
    # print('bkmer', bkmer, ' is in decset and ', 'akmer', akmer, ' is not')
    return -1


if __name__ == '__main__':
    main()
