import random

import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
import re


k = int(sys.argv[1])
L = int(sys.argv[2])


def main():
    SIZE = 1000
    s = ''
    if sys.argv[3] == 'r':
        s = np.array2string(np.random.randint(low=0, high=4, size=SIZE, dtype='i1'), separator="", suffix="",
                            prefix="")[1:-1].replace('\n', '').replace(' ', '')
    if sys.argv[3] == 'i':
        with open(sys.argv[4], 'r') as f:
            s = re.sub(r"[^ACGT]", "", f.read())
            transtable = s.maketrans('ACGT', '0123')
            s = s.translate(transtable)
            SIZE = len(s)
    # print(s)
    preds = np.zeros(pow(4, k))
    with open("../gen/xpreds/preds_" + str(k) + "_" + str(L) + ".txt", "r") as file:
        lines = file.readlines()
        for line in lines:
            transtable = line.split('\t')[0].maketrans("ACGT", "0123")
            quatstringline = line.split('\t')[0].translate(transtable)
            flot = float(line.split('\t')[1][:-1])
            preds[toDeci(quatstringline)] = flot

    with open("docksdensitysets/decyc" + str(k) + ".txt", "r") as decycfile:
        lines = decycfile.readlines()
        # print("decycling update now")
        for line in lines:
            transtable = line.replace('\n', '').maketrans("ACGT", "0123")
            quatstringline = line.translate(transtable)
            deciline = toDeci(quatstringline)
            preds[deciline] = preds[deciline] + random.uniform(0, 1)
    # print("finished decycling")
    # print("finished sorting")
    n = 1
    lidx = 0
    currmzeridx = findNewOne(s, lidx, preds)
    currmzer = s[currmzeridx: currmzeridx + k]
    lidx += 1
    while lidx < SIZE - L + 1:
        # print("lidx is ", lidx)
        if lidx - 1 == currmzeridx:
            currmzeridx = findNewOne(s, lidx, preds)
            currmzer = s[currmzeridx: currmzeridx + k]
            n += 1
        else:
            newkmeridx = lidx + L - k
            newkmer = s[newkmeridx: newkmeridx + k]
            if preds[toDeci(currmzer)] < preds[toDeci(newkmer)]:
                currmzeridx = newkmeridx
                currmzer = newkmer
                n += 1
        lidx += 1
    #print('there are ', n, ' minimizers')
    print(n / (SIZE - k + 1))


def findNewOne(s, idx, preds):
    lmer = s[idx: idx + L]
    bestKmerSoFar = 0
    bestValueSoFar = preds[toDeci(lmer[:k])]
    i = 1
    while i < L - k + 1:
        currKmer = lmer[i: i + k]
        currValue = preds[toDeci(currKmer)]
        if currValue > bestValueSoFar:
            bestValueSoFar = currValue
            bestKmerSoFar = i
        i += 1
    return idx + bestKmerSoFar


def val(c):
    if '0' <= c <= '9':
        return ord(c) - ord('0')
    else:
        return ord(c) - ord('A') + 10


def toDeci(ns):
    power = 1  # Initialize power of base
    num = 0  # Initialize result
    # print("ns is ", ns, " of length ", len(ns))
    for i in range(k - 1, -1, -1):
        # if val(ns[i]) >= 4:
        #     print('Invalid Number')
        #     return -1
        # print("i is ", i)
        num += val(ns[i]) * power
        power = power * 4
    # print("num is ", num)
    return num


if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
