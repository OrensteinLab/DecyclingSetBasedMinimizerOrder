from sys import argv
from re import sub
from randomDensityREP import findNewOne

k = int(argv[1])
L = int(argv[2])


def main():
    with open(argv[3], 'r') as f:
        s = sub(r"[^ACGTacgt]", "", f.read())
        transtable = s.maketrans('ACGTacgt', '01230123')
        s = s.translate(transtable)
    size = len(s)
    if size == 0:
        print('size 0')
        return -1
    currmzeridx = findNewOne(s, 0)
    currmzer = s[currmzeridx: currmzeridx + k]
    n = 1
    lidx = 1
    while lidx < size - L + 1:
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
    dens = n / (size - k + 1)
    print(dens)


if __name__ == '__main__':
    main()