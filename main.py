from functions import *  # import everything from functions.py

# read in data from cnt and mes files
CNT = readfile("coordinates.cnt", 1)
MES = readfile("measurements.mes", 1)

# build unknowns vector x
# x consists of (X,Y) of each unknown point
x = buildx(CNT)
# choose a-priori variance factor
sigma0 = 1
# build weight matrix P
P = buildP(MES, sigma0)

# Main loop start ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some settings for the loop. These could be put in a separate settings file
deltasum = 100
threshold = 0.0000000001
maxit = 100
count = 1
# loop until the sum of delta is less than some threshold
while deltasum>threshold:
    print("Iteration: " + str(count))  # print iteration number for the user
    # also break if maximum iterations are reached
    if count >= maxit:
        print("Maximum number of iterations reached")
        print("Break")
        break

    # build the design matrix A and the misclosure matrix w
    A, w = buildAw(CNT, MES, x)

    # calculate solution
    At = A.transpose()  # transpose design matrix
    u = At @ P @ w
    N = At @ P @ A
    Cx = np.linalg.inv(N)
    delta = -Cx @ u

    #x = x + delta
    deltasum = 0
    # for each element of delta
    for i in range(0, len(delta)):
        x[i] = x[i] + delta[i]  # add to x
        deltasum = deltasum + abs(delta[i])  # sum absolute value to deltasum

    print("deltasum = " + str(deltasum[0]) + "\n\n")
    count = count + 1  # iterate count

print('done!')