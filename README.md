# Decycling-set-based minimizer orders

This repository contains code for the manuscript [Efficient minimizer orders for large values of *k* using minimum decycling sets](https://www.biorxiv.org/content/10.1101/2022.10.18.512682)

#### Code
c++ source code for Decycling-set-based minimizers is in `src`, provided as header files:

`minimzers.h`: decycling set based minimizers (up to *k*=63)

`big_minimzers.h`: decycling set based minimizers for *k* longer than 63. This uses the [GMP bignum library](https://gmplib.org/) for larger integer representations. You must install the GMP library and link it when compiling your code to use these minimizers.

`test_minimizers.cpp`: an example of how to use minimizers in your code. It can be used to compute the density for an input sequence or random sequence.

#### Compilation 
`g++ -Wall -O3 -o test_minimizer test_minimizers.cpp -std=c++2a -lgmp -lgmpxx` Note that `-lgmp -lgmpxx` are **only** necessary for big_minimizers

#### Usage
Command line options for test_minimizer:
```
Usage: ./test_minimizer -k <k> -L <L> -r/--reps -l/--seqlen -m/--method -seq/--seqfile -o/--outfile [-set/--setfile]
        k:       k-mer length
        L:       window length
        reps:    # of randomized repeats
        seqlen:  Length of random sequence to generate.
                 If input sequence is read from file and length is non-zero, this length will be sampled from the input.
        method:  Order to use. Options: 'random', 'docks', 'pasha', 'uhs', 'decycling', 'double', 'set'
        seqfile: The sequence is read from this file, or if it doesn't exist, the random sequence will be saved to this file
        outfile: File to append output statistics to.
        setfile: The file with the set to use. Required when method is 'set', 'uhs', 'pasha', or 'docks'
```

#### Other files
`results` contains tables of raw results reported in the paper

`scripts/all_plots.py` contains plotting code to generate the paper figures from the results tables


