# DecyclingSetBasedMinimizerOrder

DensityKL.py command structure info:

python3 DensityKL.py [order type flag] [kmin] [kmax] [lmin] [lmax] [repetitions] [input sequence - 'i'/random sequence - 'r'] [if 'i' then pathfile]

order types flags:
'-l' = lexicographical
'-r' = random
'-d' = DOCKS-based (random internal order)
'-dl' = DOCKS-basd (lexicographic internal order)
'-p' = prediction-based
'-c' = decycling/random-based
'-u' = UHS-based, decycling preferred
'-s' = decycling/symmetric
'-cm' =  modified decycling
'-sm' = modified decycling/symmetric

command example:
'python3 DensityKL.py -d 8 12 20 100 5 i sequence.txt'

-------------------------------------------------------------

plotting.py command structure info:

python3 plotting.py [graph dependence flag (-k or -l)] [constant value of k or l (if previous flag is '-k' than it's value of l and vice versa)] [flags indicating orders types to show in graph] [input sequence flag == 'is'] [name of sequence (only if previous flag is on)] [directory path of the source tables]

order type flags are the same as shown above but without the '-'

command example:
'python3 plotting.py -k 100 d c u p is ch1HomoSap ./tables/ch1HomoSap/k_from50to200'

this will generate a graph where y is density x is k value that shows plots of order types: docks, decycling, uhs based, predictions. of the input sequence ch1HomoSap, from the tables in './tables/ch1HomoSap/k_from50to200'.
