from functions import * # import everything from functions.py
from pandas import * # used for debugging only

# read in data from cnt and mes files
CNT = readfile("coordinates.cnt", 1)
MES = readfile("measurements.mes", 1)

# build x vector
# x consists of (X,Y) of each unknown point
x = buildx(CNT)
print(x)