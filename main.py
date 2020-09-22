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

    # Main loop START ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # some settings for the loop. These could be put in a separate settings file
    deltasum = 100  # initialize deltasum to some large value
    threshold = 0.0000000001  # threshold for minimum deltasum
    # this threshold is way smaller than it needs to be, but since the function converges very quickly this doesn't noticeably impact performance
    maxit = 100  # maximum iterations before breaking the loop
    count = 1  # initialize iteration counter
    # loop until the absolute sum of the elements in delta is less than some small threshold
    while deltasum > threshold:
        if not suppress_print:
            print("Iteration: " + str(count))  # print iteration number for the user
        # also break if maximum number of iterations is reached
        if count >= maxit:
            print("Maximum number of iterations reached")
            print("Break")
            break

        # build the design matrix A and the misclosure vector w
        A, w = buildAw(CNT, MES, x)

        # calculate solution
        At = A.transpose()  # transpose design matrix
        u = At @ P @ w
        N = At @ P @ A
        Cx = np.linalg.inv(N)
        delta = -Cx @ u

        deltasum = 0
        # for each element of delta
        for i in range(0, len(delta)):
            x[i] = x[i] + delta[i]  # add to x
            deltasum = deltasum + abs(delta[i])  # sum absolute value to deltasum

        if not suppress_print:
            print("deltasum = " + str(deltasum[0]) + "\n")
        count = count + 1  # iterate count
    # Main loop END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate residuals
    rhat = A @ delta + w
    # Adjusted observations
    lhat = np.zeros([len(MES), 1])
    for i in range(0, len(MES)):
        lhat[i] = MES[i][3] + rhat[i]
    # Calculate posteriori variance factor and unit variance factor
    sigma0hat = (rhat.transpose() @ P @ rhat)/(len(MES) - len(x))
    sigma0hat = sigma0hat[0][0]  # convert to float from ndarray
    unitvar = sigma0hat/sigma0
    # Calculate the v-c matrix of unknowns Cxhat
    Cxhat = sigma0hat * Cx
    # Calculate the v-c matrix of measurements Clhat
    Clhat = A @ Cxhat @ At
    # Calculating the v-c matrix of residuals Crhat
    Crhat = sigma0hat * np.linalg.inv(P) - Clhat
    # Calculate error ellipse
    # for each 2x2 sub-matrix in Cxhat
    for i in range(0, len(Cxhat), 2):
        # extract sub matrix
        subM = np.zeros([2, 2])
        subM[0][0] = Cxhat[i][i]
        subM[0][1] = Cxhat[i][i+1]
        subM[1][1] = Cxhat[i+1][i+1]
        subM[1][0] = Cxhat[i+1][i]



    if not suppress_print:
        print('posteriori variance factor: ' + str(sigma0hat)
              + '\nunit variance factor: ' + str(unitvar)
              )
        print('done!')


# run main function
if __name__ == '__main__':
    main("coordinates.cnt", "measurements.mes")
