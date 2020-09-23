import timeit
import functools
from main import *

n = 10000
total = timeit.timeit(functools.partial(main, "coordinates.cnt", "measurements.mes", suppress_print=True, plot=False), setup='gc.enable()', number=n)
print("Average Time: " + str(total/n) + " seconds")
