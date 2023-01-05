# import numpy
# import pandas as pd

# import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from sys import argv

NUMOFORDERTYPES = 10
# NUMOFLVALUES = 10
# NUMOFKVALUES = 6
SEQSIZE=1000000

def main():

    mode = 0
    expected = True
    if argv[1] == '-k':
        mode = 'k'
    elif argv[1] == '-l':
        mode = 'l'
    else:
        print('invalid mode')
        exit(-1)

    nofplots = 0
    whichplots = np.full((NUMOFORDERTYPES), False)

    plotdic = {0: 'lex', 1: 'random', 2: 'docks', 3: 'randec',
                        4: 'uhs', 5: 'preds', 6: 'dockslex', 7: 'decsym', 8: 'decsymmod', 9: 'decmod'}

    colors = ['red', 'black', 'blue', 'darkgreen', 'rebeccapurple', 'lime', 'darkorange', 'deeppink', 'sienna', 'lightcoral']
    # markers = ['p', 's', 'o', 'D', '|', '^', 'x', '+', '*', ',']

    for x in argv:
        if x == 'l':
            whichplots[0] = True
            nofplots += 1
        if x == 'r':
            whichplots[1] = True
            nofplots += 1
        if x == 'd':
            whichplots[2] = True
            nofplots += 1
        if x == 'c':
            whichplots[3] = True
            nofplots += 1
        if x == 'u':
            whichplots[4] = True
            nofplots += 1
        if x == 'p':
            whichplots[5] = True
            nofplots += 1
        if x == 'dl':
            whichplots[6] = True
            nofplots += 1
        if x == 's':
            whichplots[7] = True
            nofplots += 1
        if x == 'sm':
            whichplots[8] = True
            nofplots += 1
        if x == 'cm':
            whichplots[9] = True
            nofplots += 1
        if x == 'is':
            expected = False

        if x == 'all':
            whichplots = np.where(whichplots > -1, True, False)
            break

    marks = np.nonzero(whichplots)[0]

    plotlist, stdlist = getTables(marks, plotdic)
    names = np.array(
        ['Lexicographic', 'Random', 'UHS', 'Decycling', 'Decycling-UHS',
         'Predictions-based', 'DOCKS-based(lex. internal order)', 'Double Decycling', 'Modified Double Decycling', 'Modified Decycling'])
    legend = []

    korl = argv[2]
    min = 100
    max = 0
    if mode == 'k':
        l = int(argv[2])
    else:
        k = int(argv[2])
        for t in range(0, NUMOFORDERTYPES):
            if plotlist[t] is not None:
                plotlist[t] = plotlist[t].transpose()
                stdlist[t] = stdlist[t].transpose()
    i = 0
    for x in marks:
        xvals = plotlist[x][0][1:]
        if mode == 'k':
            print(f'plotlist[x][:,0] is {plotlist[x][:,0]}')
            row = np.where(plotlist[x][:, 0] == int(korl))[0][0]
        else:
            row = np.where(plotlist[x][:, 0] == int(korl))[0][0]
        ydenvals = plotlist[x][row][1:]
        ystdvals = stdlist[x][row][1:]
        currmin = xvals[0]
        currmax = xvals[-1]
        if currmin < min:
            min = currmin
        if currmax > max:
            max = currmax

        if mode == 'k':
           factors = (l-xvals+2)*(SEQSIZE-xvals+1)/(SEQSIZE-l+1)
        else:
           factors = (xvals-k+2)*(SEQSIZE-k+1)/(SEQSIZE-xvals+1)

        print(factors)
        plt.plot(xvals, ydenvals*factors, color=colors[i])#, marker=markers[x])
        # print('xvals: ', xvals)
        # print('ystdvals: ', ystdvals)
        # print('ydenvals: ', ydenvals)
        plt.errorbar(xvals, ydenvals*factors, yerr=ystdvals*factors, ls='None', color=colors[i])
        i += 1
        legend.append(names[x])
    if mode == 'k':
        xtics = np.arange(min, max+1, 1)
    else:
        xtics = np.arange(min, max + 1, 10)
    title = "Density Factor "
    if expected:
        title = 'Expected ' + title
    else:
        title += 'on CHM13X '
    if mode == 'k':
        title += 'for L=' + str(korl)
        plt.xlabel("K")
    else:
        title += 'for K=' + str(korl)
        plt.xlabel("L")
     #   plt.yscale("log")
    plt.title(title)
    plt.xticks(xtics)
    # plt.xlabel("K")
    if expected:
        plt.ylabel("Expected Density")
    else:
        plt.ylabel("Particular Density")
    plt.legend(legend)
    # plt.show()
    plt.savefig("./comparison graphs/" + title + ".pdf")


def getTables(marks, dict):
    plotlist = [None, None, None, None, None, None, None, None, None, None]
    stdlist = [None, None, None, None, None, None, None, None, None, None]
    for x in marks:
        plotlist[x] = genfromtxt(argv[-1] + dict[x] + 'Density.csv', delimiter=' ')
        stdlist[x] = genfromtxt(argv[-1] + dict[x] + 'STD.csv', delimiter=' ')

    return plotlist, stdlist


if __name__ == '__main__':
    main()
