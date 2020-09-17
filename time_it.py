import timeit
import functools
from main import *

n = 10000
total = timeit.timeit(functools.partial(main, True), setup='gc.enable()', number=n)
print("Average Time: " + str(total/n) + " seconds")
