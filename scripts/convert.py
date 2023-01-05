import numpy as np
import pandas as pd
import itertools
import os,sys,glob

indir = sys.argv[1]

basedict={"A": 0, "C": 1, "G": 2, "T": 3}

decfiles = glob.glob(indir+'/decyc11.txt')

for decfile in decfiles:
  k = decfile.split('decyc')[1].split(".")[0]
  setfiles = glob.glob(indir+'/*hit{}_*.txt'.format(k))
  
  decint_file = indir+'/decyc{}.int'.format(k)
  with open(decfile) as f, open(decint_file,'w') as o:
      for line in f:
          kmer = line.strip("\n")
          totval = 0
          for pos,char in enumerate(kmer):
             val = basedict[char] * pow(4, (int(k) - pos - 1))
             totval += val
          o.write(str(totval) + "\n")
  for setfile in setfiles:
     L = setfile.split('_')[-1].split('.txt')[0]
     set_intfile = indir + "/hit_{}_{}.int".format(k,L)
     with open(decfile) as f, open(set_intfile,'w') as o:
       for line in f:
          kmer = line.strip("\n")
          totval = 0
          for pos,char in enumerate(kmer):
             val = basedict[char] * pow(4, (int(k) - pos - 1))
             totval += val
          o.write(str(totval) + "\n")
     with open(setfile) as f, open(set_intfile,'a') as o:
       for line in f:
          kmer = line.strip("\n")
          totval = 0
          for pos,char in enumerate(kmer):
             val = basedict[char] * pow(4, (int(k) - pos - 1))
             totval += val
          o.write(str(totval) + "\n")
      
  

