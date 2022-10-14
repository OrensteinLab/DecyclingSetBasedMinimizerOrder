import numpy as np
import sys
import os
import os.path
import subprocess
import re
from subprocess import check_output
from time import time

np.set_printoptions(threshold=sys.maxsize)
kmin = int(sys.argv[2])
kmax = int(sys.argv[3])
lmin = int(sys.argv[4])
lmax = int(sys.argv[5])
REP = int(sys.argv[6])


def main():
    prog = ''
    filename = ''
    seqpath = ''
    seqname = ''
    rows = int((lmax - lmin + 20) / 10)
    cols = kmax - kmin + 2
    densitytable = np.zeros([rows, cols])
    stdtable = np.zeros([rows, cols])
    mode = sys.argv[7]
    for llabel in range(1, rows):
        densitytable[llabel][0] = int(lmin + 10 * (llabel - 1))
        stdtable[llabel][0] = int(lmin + 10 * (llabel - 1))
    for klabel in range(1, cols):
        densitytable[0][klabel] = int(kmin + klabel - 1)
        stdtable[0][klabel] = int(kmin + klabel - 1)

    if sys.argv[1] == '-l':
        prog = 'lexDensityREP'
        filename = 'lex'
    elif sys.argv[1] == '-r':
        prog = 'randomDensityREP'
        filename = 'random'
    elif sys.argv[1] == '-d':
        prog = 'docksDensityREP'
        filename = 'docks'
    elif sys.argv[1] == '-dl':
        prog = 'dockslexDensity'
        filename = 'dockslex'
    elif sys.argv[1] == '-p':
        prog = 'predsDensity'
        filename = 'pred'
    elif sys.argv[1] == '-c':
        prog = 'randecDensityNB'
        filename = 'randec'
    elif sys.argv[1] == '-u':
        prog = 'uhsDensityMP'
        filename = 'uhs'
    elif sys.argv[1] == '-s':
        prog = 'decsymDensityNB'
        filename = 'decsym'
    elif sys.argv[1] == '-cm':
        prog = 'randecDensityMOD'
        filename = 'randecMOD'
    elif sys.argv[1] == '-sm':
        prog = 'decsymDensityMOD'
        filename = 'decsymMOD'

    if mode == 'i':
        seqpath = sys.argv[8]
        seqname = seqpath.split('/')[-1].split('.')[0]

    for k in range(kmin, kmax + 1):
        for l in range(lmin, lmax + 10, 10):

            command = 'python3 ' + prog + '.py ' + str(k) + ' ' + str(l)
            command += ' ' + str(REP)
            if mode == 'i':
                command = command + ' i ' + seqpath
            elif mode == 'r':
                command = command + ' r'
            densityrunstart = time()
            output = str(check_output(command, shell=True))
            densityrunend = time()
            print(f'output is {output} for ({k}, {l}) - took {densityrunend - densityrunstart} seconds')
            mean = output.split(' ')[0]
            std = output.split(' ')[1]
            mean = re.sub(r"[^\d|.]", "", mean)
            std = re.sub(r"(n|\\|\')", "", std)
            print(f'mean:{mean}, std:{std}')
            colidx = k - kmin + 1
            rowidx = int((l - lmin + 10)/10)
            densitytable[rowidx][colidx] = float("{0:.6f}".format(float(mean)))
            stdtable[rowidx][colidx] = float(std)
            # print("density for ", k, "/", l, " is ", mean, " with std ", std)
    print(f'densitytable: {densitytable}')
    dirname = './tables' + '/' + seqname
    # print('dirname is ' + dirname)
    if not os.path.isdir(dirname):
        print('no dir named ' + dirname + ' so creating one')
        p = subprocess.Popen('mkdir ' + dirname, shell=True)
        p.wait()

    # print('saving to ' + dirname + '/ bc is exists? ' + str(os.path.isdir(dirname)))
    np.savetxt(dirname + '/' + prog + 'K' + str(kmin) + 'to' + str(kmax) + 'L' + str(lmin) + 'to' + str(lmax) + '.csv', densitytable, delimiter=' ')
    np.savetxt(dirname + '/' + filename + 'STD' + 'K' + str(kmin) + 'to' + str(kmax) + 'L' + str(lmin) + 'to' + str(lmax) + '.csv', stdtable, delimiter=' ')


# def density(cmd):
#     outputs = np.empty([REP])
#     for i in range(REP):
#         output = check_output(cmd, shell=True)
#         # stream = os.popen(cmd)
#         # output = stream.read()
#         print('output is', output)
#         outputs[i] = float(output)
#
#     output = np.empty([2])
#     output[0] = np.mean(outputs)
#     output[1] = np.std(outputs)
#     # print(" output is ", output)
#     return output




if __name__ == '__main__':
    main()
