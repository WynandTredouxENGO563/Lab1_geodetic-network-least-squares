import timeit
import functools
from main import *
import sys

def benchmark(n=10000):
    total = timeit.timeit(functools.partial(main, "coordinates.cnt", "measurements.mes", suppress_print=True, plot=False), setup='gc.enable()', number=n)
    print("Total Time: " + str(total) + " seconds\nAverage Time: " + str(total/n) + " seconds")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        benchmark()
    else:
        benchmark(int(sys.argv[1]))
