from functions import *  # import everything from functions.py


# main function
# CNTfile: Filename of coordinates file
# MesFile: Filename of measurements file
# CNTheader: number of header lines in coordinates file (default is 1)
# MESheader: number of header lines in measurements file (default is 1)
# suppress_print: Set to True to suppress outputs to the console
def main(CNTfile, MesFile, CNTheader=1, Mesheader=1, suppress_print=False):
    # read in data from cnt and mes files
    CNT = readfile(CNTfile, CNTheader)
    MES = readfile(MesFile, Mesheader)

    # build unknowns vector x
    # x consists of (X,Y) of each unknown point
    x = buildx(CNT)
    # choose a-priori variance factor
    sigma0 = 1
    # build weight matrix P
    P = buildP(MES, sigma0)

    # Main loop start ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # some settings for the loop. These could be put in a separate settings file
    deltasum = 100  # initialize deltasum to some large value
    threshold = 0.0000000001  # threshold for minimum deltasum
    maxit = 100  # maximum iterations before breaking the loop
    count = 1 # initialize iteration counter
    # loop until the sum of delta is less than some threshold
    while deltasum > threshold:
        if not suppress_print:
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

        if not suppress_print:
            print("deltasum = " + str(deltasum[0]) + "\n\n")
        count = count + 1  # iterate count

    if not suppress_print:
        print('done!')


# run main function
if __name__ == '__main__':
    main("coordinates.cnt", "measurements.mes")
