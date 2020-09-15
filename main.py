from functions import * # import everything from functions.py
from pandas import * # used for debugging only

CNT = readfile("coordinates.cnt", 1)
print(DataFrame(CNT))
print()
MES = readfile("measurements.mes", 1)
print(DataFrame(MES))