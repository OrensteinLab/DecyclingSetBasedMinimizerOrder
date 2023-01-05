from numpy import zeros, dot
from itertools import starmap, repeat
from ctypes import c_uint32


def bitstoArray(k, set):
    iterations = max(int(pow(4, k) / 32), 1)
    arraybits = zeros(iterations, dtype=c_uint32)
    list(starmap(buildArray, zip(range(0, iterations), repeat(set), repeat(arraybits))))
    return arraybits


def buildArray(iter, barray, arraybits):
    # print(f'barray is {barray}')
    # print(f'iter is {iter}')
    bitidx = iter << 5
    section = barray[bitidx: bitidx + 32]
    # print(f'section is {section}')
    st = section.to01()
    # print(f'st is {st}')
    toadd = int(st, 2)
    # print(f'toadd is {toadd}')
    arraybits[iter] = toadd


def checkmembership(kmer, array):
    idx = toDeci(kmer)
    arid = int(idx / 32)
    internalid = idx % 32
    setnumber = format(array[arid], '#034b')[2:]
    return int(setnumber[internalid])


def checkIndec(kmer, coeffs):  # weights,arange):
    k = len(kmer)
    intkmer = [int(x) for x in kmer]
    # print(f'kmer:{intkmer}||coeffs:{coeffs}')
    s = dot(intkmer[1:], coeffs)
    #   s = sum(weights[arange,kmer[1:]])
    if s < -0.00000001:
        return 0
    elif s > 0.00000001:
        shiftsum = dot(intkmer[:-1], coeffs)
        #       shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum > 0.00000001:
            return 0
        else:
            return 1
    else:
        shiftsum = dot(intkmer[:-1], coeffs)
        #       shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum > 0.00000001 or shiftsum < -0.00000001:
            return 0
        for i in range(1, k):
            didbreak = False
            j = 0
            for j in range(k - i):
                if intkmer[j] > intkmer[(i + j)]:
                    return 0
                if intkmer[j] < intkmer[(i + j)]:
                    didbreak = True
                    break
            if j == k - 2: return 1 # homopolymer
            if not didbreak:  # got here by finishing inner loop
                if k % i == 0:
                    return 1
                else:
                    return 0
        return 1


def checkInsym(kmer, coeffs):
    k = len(kmer)
    intkmer = [int(x) for x in kmer]
    s = dot(intkmer[1:], coeffs)
    # s = sum(weights[arange,kmer[1:]])
    if s > 0.000001:
        return False
    elif s < -0.000001:
        shiftsum = dot(intkmer[:-1], coeffs)
        #    shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum < -0.000001:
            return False
        else:
            return True
    else:
        shiftsum = dot(intkmer[:-1], coeffs)
        #    shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum > 0.000001 or shiftsum < -0.000001:
            return False
        for i in range(1, k):
            didbreak = False
            j = 0
            for j in range(k - i):
                if intkmer[j] > intkmer[(i + j)]:
                    return False
                if intkmer[j] < intkmer[(i + j)]:
                    didbreak = True
                    break
            if j == k - 2: return True  # homopolymer
            if not didbreak:  # got here by finishing inner loop
                if k % i == 0:
                    return True
                else:
                    return False
        return True


def checkInsymmod(kmer, coeffs):
    intkmer = [int(x) for x in kmer]
    s = dot(intkmer[1:], coeffs)
    #   s = sum(weights[arange,kmer[1:]])
    if s > -0.000001:
        return False
    else:
        shiftsum = dot(intkmer[:-1], coeffs)
        # shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum < -0.000001:
            return False
        else:
            return True


def checkIndecmod(kmer, coeffs):
    intkmer = [int(x) for x in kmer]
    s = dot(intkmer[1:], coeffs)
    #  s = sum(weights[arange,kmer[1:]])
    if s < 0.000001:
        return False
    else:
        shiftsum = dot(intkmer[:-1], coeffs)
        #      shiftsum = sum(weights[arange,kmer[:-1]])
        if shiftsum > 0.000001:
            return False
        else:
            return True


def val(c):
    if '0' <= c <= '9':
        return ord(c) - ord('0')
    else:
        return ord(c) - ord('A') + 10


def toDeci(ns):
    k = len(ns)
    power = 1
    num = 0
    # print("ns is: .", ns, ".")
    for i in range(k - 1, -1, -1):
        # if val(ns[i]) >= 4:
        #     print('Invalid Number')
        #     return -1
        # print("i is ", i)
        num += val(ns[i]) * power
        power = power * 4
    # print('num is ', num)
    return num


def uhsfindNewOne(s, idx, k, L, compareFunc, coeffs, uhsset):
    lmer = s[idx: idx + L]
    bestKmerSoFaridx = 0
    bestkmersofar = lmer[:k]
    i = 1
    while i < L - k + 1:
        currKmer = lmer[i: i + k]
        if compareFunc(currKmer, bestkmersofar, uhsset, coeffs) > 0:
            bestkmersofar = currKmer
            bestKmerSoFaridx = i
        # else:
        #     print('findnewone: not chosen')
        i += 1
    return idx + bestKmerSoFaridx


def decfindNewOne(s, idx, k, L, compare, coeffs):
    lmer = s[idx: idx + L]
    bestKmerSoFaridx = 0
    bestkmersofar = lmer[:k]
    i = 1
    while i < L - k + 1:
        currKmer = lmer[i: i + k]
        if compare(currKmer, bestkmersofar, coeffs) > 0:
            bestkmersofar = currKmer
            bestKmerSoFaridx = i
        i += 1
    return idx + bestKmerSoFaridx


