from functions import *  # import everything from functions.py
from datetime import datetime


# main function
# CNTFile: Filename of coordinates file
# MESFile: Filename of measurements file
# CNTheader: number of header lines in coordinates file (default is 1)
# MESheader: number of header lines in measurements file (default is 1)
# suppress_print: Set to True to suppress outputs to the console
# plot: Set to False to suppress plotting
# sigma0: a-priori variance factor. Chosen as 1 by default
# maxit: maximum iterations before breaking the loop
def main(CNTFile='', MESFile='', CNTheader=1, Mesheader=1, suppress_print=False, plot=True, sigma0=1, maxit=100):
    # start timer
    time0 = datetime.now()
    # if no CNT filename was provided
    if CNTFile == '':
        # automatically find CNT file in directory
        CNTFile = FindFiles('cnt')
    # if no MES filename was provided
    if MESFile == '':
        # automatically find CNT file in directory
        MESFile = FindFiles('mes')

    # read in data from cnt and mes files
    CNT = readfile(CNTFile, CNTheader)
    MES = readfile(MESFile, Mesheader)

    # build unknowns vector x
    # x consists of (X,Y) of each unknown point
    x = buildx(CNT)
    # build weight matrix P
    P = buildP(MES, sigma0)

    # Main loop START ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # some settings for the loop. These could be put in a separate settings file
    deltasum = 100  # initialize deltasum to some large value
    # threshold for minimum deltasum should be 1/2 of the smallest measurement standard deviation
    threshold = 0.5*min([i[4] for i in MES])
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
    time1 = datetime.now()
    timetaken = time1 - time0

    # plot points using matplotlib
    if plot:
        main_ax = plotCNT(CNT, x)

    # Open output file
    out = open("output.out", "w")
    divider = '\n\n' + '*'*100 + '\n'
    line = '\n' + '-'*50
    # count number of angle/dist measurements
    num_angle = sum([i[2] == 'Angle' for i in MES])
    num_dist = sum([i[2] == 'Dist' for i in MES])
    out.write(
"""Geodetic Network Least Squares Adjustment
Wynand Tredoux -- September 2020\n
Execution Date:\t""" + str(time0) + """
Execution Time:\t""" + str(timetaken.total_seconds()) + """ seconds
Itterations:\t""" + str(count-1) + """
Threshold:\t\t""" + str(threshold) + """
sigma0:\t\t\t""" + str(sigma0) + divider + """
Observations/Unknowns Summery\n
Angle Measurements:\t\t\t""" + str(num_angle) + """
Distance Measurements:\t\t""" + str(num_dist) + line + """
Total Measurements:\t\t\t""" + str(len(MES)) + """
Total unknowns:\t\t\t\t""" + str(len(x)) + line + """
Total Degrees of Freedom:\t""" + str(len(MES) - len(x)) + """
\n"""
    )

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

    # Write estimates unknowns to output file
    out.write(
"""Posteriori Variance Factor:\t""" + str(sigma0hat) + """
Unit Variance Factor:\t\t""" + str(unitvar) + divider + """
Estimated Unknowns\n
Point Name\tx\ty\tdiffx\tdiffy"""
    )
    # get old unknown points from CNT
    # compare with estimates coordinates
    max_name_len = max([len(i[0]) for i in CNT]) + 4  # maximum name length in CNT plus 4 buffer spaces
    count = 0
    unknown_points = []
    for i in range(0,len(CNT)): # for all points in CNT
        P = Point(CNT[i])
        if P.isUnknown(): # if point is unknown
            unknown_points.append(P)
            # get estimates x and y from x
            xest = x[count*2]
            yest = x[count*2 + 1]
            # calc estimated(x,y) - initial(x,y)
            diffx = xest - P.x
            diffy = yest - P.y
            max_num_len = max(numlen([P.x, P.y, diffx, diffy])) # find maximum length of numbers not including decimals
            decimals = 4  # number of decimal places to keep
            # print to output file
            outstr = '\n{:' + str(max_name_len) + '}'
            # spacing should be max_num_len + decimals + 2 for decimal point and signs + 2 for buffer
            tmp = '{:' + str(max_num_len + decimals + 4) + '.' + str(decimals) + 'f}'
            outstr = outstr + tmp*4
            outstr = outstr.format(P.name, xest, yest, diffx, diffy)
            out.write(outstr)
            count = count + 1

    # write rhat to file
    out.write(divider + "\nVector of Residuals\n")
    np.savetxt(out, rhat, fmt='%+1.6E')
    # Write lhat to file
    out.write(divider[1:] + "\nVector of Corrected Measurements\n")
    np.savetxt(out, lhat, fmt='%12.6f')
    # Write Cxhat to file
    out.write(divider[1:] + "\nVariance Covariance Matrix of Unknowns\n")
    np.savetxt(out, Cxhat, fmt='%+1.4E', delimiter='\t')
    # Write Clhat to file
    out.write(divider[1:] + "\nVariance Covariance Matrix of Corrected Measurements\n")
    np.savetxt(out, Clhat, fmt='%+1.4E', delimiter='\t')
    # Write Crhat to file
    out.write(divider[1:] + "\nVariance Covariance Matrix of Residuals\n")
    np.savetxt(out, Crhat, fmt='%+1.4E', delimiter='\t')

    # Calculate error ellipse
    out.write(divider[1:] + """\nError Ellipse\n
Name\tSemi-Major axis\tSemi-Minor axis\ttheta"""
              )
    # for each 2x2 sub-matrix in Cxhat
    for i in range(0, len(Cxhat), 2):
        # extract sub matrix
        subM = np.zeros([2, 2])
        subM[0][0] = Cxhat[i][i]
        subM[0][1] = Cxhat[i][i+1]
        subM[1][1] = Cxhat[i+1][i+1]
        subM[1][0] = Cxhat[i+1][i]
        # calc eigen values
        evals = np.linalg.eigvals(subM)
        lmax = evals.max()
        lmin = evals.min()
        # calc semi-minor/major axis'
        semi_minor = math.sqrt(lmin)
        semi_major = math.sqrt(lmax)
        # calc orientation of semi-major axis
        theta = 0.5 * math.atan2(-2*subM[0][1], subM[1][1] - subM[0][0])  # theta is the angle counter-clockwise from the x axis to the semi-major axis

        # write to file
        outstr = '\n{:' + str(max_name_len) + '}'
        max_num_len = max(numlen([theta, lmax, lmin]))  # find maximum length of numbers not including decimals
        decimals = 8  # number of decimal places to keep
        # spacing should be max_num_len + decimals + 2 for decimal point and signs + 3 for buffer
        tmp = '{:' + str(max_num_len + decimals + 5) + '.' + str(decimals) + 'f}'
        # spacing for E should be decimals + 1 for E + 1 for sign + 2 for exponent + 2 for number and decimal + 4 for buffer
        tmp2 = '{:' + str(decimals + 10) + '.' + str(decimals) + 'E}'
        outstr = outstr + tmp2*2 + tmp*1
        outstr = outstr.format(unknown_points[int(i/2)].name, semi_major, semi_minor, theta)
        out.write(outstr)

        # draw error ellipses on figure f1
        if plot:
            drawEE(main_ax, unknown_points[int(i/2)].x, unknown_points[int(i/2)].y, semi_major, semi_minor, theta, scale=10000, name=unknown_points[int(i/2)].name)
    # add legend to f1
    if plot:
        plt.figure(main_ax.figure.number)  # make main_ax the active window
        legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#fc4c4c', label='Unknown'),
                           Line2D([0], [0], marker='o', color='w', markerfacecolor='#03c2fc', label='Known'),
                           pat.Ellipse(xy=(0, 0), width=0, height=0,
                                       angle=0, facecolor='none', edgecolor='red', label='Error Ellipse (exagerated)')]
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.1, 1), loc='upper right')

    if not suppress_print:
        print('posteriori variance factor: ' + str(sigma0hat)
              + '\nunit variance factor: ' + str(unitvar)
              )
        print('done!')
    out.close()
    # save figures to pdf
    if plot:
        SaveFigs("Figures.pdf")
    # show the figure (this also pauses the program which is why it is last)
    if plot:
        plt.show()
        # close all
        plt.close('all')



# run main function
if __name__ == '__main__':
    main()
