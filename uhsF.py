from sys import argv
from re import sub
from math import sin
from multiprocessing import Array
from bitarray import bitarray
import functions
from uhsDensityMP import compare
from ctypes import c_uint32

k = int(argv[1])
L = int(argv[2])


def main():
    global coeffs
    u = 3.141592653589793 * 2.0 / k
    coeffs = [sin(i * u) for i in range(1, k)]
    with open(argv[3], 'r') as f:
        s = sub(r"[^ACGTacgt]", "", f.read())
        transtable = s.maketrans('ACGTacgt', '01230123')
        s = s.translate(transtable)
    size = len(s)
    if size == 0:
        print('size 0')
        return -1
    uhsset = bitarray('0' * pow(4, k))
    uhspath = "./docksdensitysets/PASHA" + str(k) + "_" + str(L) + ".txt"
    with open(uhspath, "r") as uhsfile:
        lines = uhsfile.readlines()
        for line in lines:
            transtable = line.maketrans("ACGT", "0123")
            translatedline = line.translate(transtable).replace('\n', '')
            idx = functions.toDeci(translatedline)
            uhsset[idx] = 1
    uhsarr = Array(c_uint32, functions.bitstoArray(k, uhsset))
    del uhsset
    currmzeridx = functions.uhsfindNewOne(s, 0, k, L, compare, coeffs, uhsarr)
    currmzer = s[currmzeridx: currmzeridx + k]
    n = 1
    lidx = 1
    while lidx < size - L + 1:
        if lidx - 1 == currmzeridx:
            currmzeridx = functions.uhsfindNewOne(s, lidx, k, L, compare, coeffs, uhsarr)
            currmzer = s[currmzeridx: currmzeridx + k]
            n += 1
        else:
            newkmeridx = lidx + L - k
            newkmer = s[newkmeridx: newkmeridx + k]
            if compare(newkmer, currmzer, uhsarr, coeffs) > 0:
                currmzeridx = newkmeridx
                currmzer = newkmer
                n += 1
        lidx += 1
    dens = n / (size - k + 1)
    print(dens)


if __name__ == '__main__':
    main()