from functions import *  # import everything from functions.py

# read in data from cnt and mes files
CNT = readfile("coordinates.cnt", 1)
MES = readfile("measurements.mes", 1)
n = len(MES)  # n = number of measurements

# build unknowns vector x
# x consists of (X,Y) of each unknown point
x = buildx(CNT)
u = len(x)  # u = number of unknowns
# choose a-priori variance factor
sigma0 = 1
# build weight matrix P
P = buildP(MES, sigma0)
# build the design matrix A and the misclosure matrix w
A, w = buildAw(CNT, MES, x, P)



print('done!')