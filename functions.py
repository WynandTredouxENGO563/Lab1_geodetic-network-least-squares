from classes import *
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import os

# Function to read in text files for this project
# returns either a string (for .txt files) or a 2D list (for .cnt or .mes files)
# Filesname: name of file including extension
# header_lines: number of lines that the header takes up (will be skipped)
def readfile(filename, header_lines):
    # acceptable file extensions
    ext = ['txt', 'mes', 'cnt']
    # get file extension
    tmp = filename.split('.')
    # check if file extension is acceptable
    if not any(x in tmp[-1] for x in ext):
        # if file is of the wrong type, raise exception
        exception_text = "File must have one of the following extensions:"
        for i in ext:
            exception_text = exception_text + "\n." + i
        raise Exception(exception_text)
    # save file extension
    filetype = tmp[-1]

    # read in file
    f = open(filename, "r")
    text = f.read()
    f.close()
    # object that will be returned
    output = []

    # if file is a regular txt file, just pass the raw text as a string
    if filetype == 'txt':
        output = text
    # if file is of one of the other types
    else:
        # split by line
        lines = text.split('\n')
        # loop through each line, ignoring headers
        linenum = header_lines  # keep track of which line we're on
        for line in lines[header_lines:]:
            linenum = linenum + 1
            # skip line if empty
            if not line:
                continue
            # split line by tab
            elements = line.split('\t')
            # parse/convert elements according to file extension
            exception_text = 'Error on line ' + str(linenum) + ' in file ' + filename  # used in Exception() if there's an error

            # for cnt files
            if filetype == 'cnt':
                # cnt files should have 4 elements per row
                if len(elements) != 4:
                    raise Exception(exception_text)
                # try to convert elements in row to proper type
                try:
                    row = [elements[0], elements[1], float(elements[2]), float(elements[3])]
                except:
                    raise Exception(exception_text)
                # save row in output
                output.append(row)

            # for mes files
            elif filetype == 'mes':
                # cnt files should have 4 elements per row
                if len(elements) != 5:
                    raise Exception(exception_text)
                # try to convert elements in row to proper type
                try:
                    # convert DMS to radians for angle measurements
                    if elements[2] == "Angle":
                        DMS = elements[3].split(' ')
                        elements[3] = (float(DMS[0]) + float(DMS[1]) / 60 + float(DMS[2]) / 3600) * math.pi/180
                        # also convert standard deviation from arcseconds to rads
                        elements[4] = float(elements[4])/3600 * math.pi/180

                    row = [elements[0], elements[1], elements[2], float(elements[3]), float(elements[4])]
                except:
                    raise Exception(exception_text)
                # save row in output
                output.append(row)

    return output


# Function to build the unknowns vector x
# returns x as a numpy array
# CNT: 2D list containing data from the coordinates.cnt file
def buildx(CNT):
    # loop through CNT and find points of type U (unknowns)
    x = np.array([])
    for row in CNT:
        type = row[1]
        if type == 'U' or type == 'u':
            x = np.append(x, [row[2], row[3]])

    return x


# Function to build the weight matrix P
# returns P as a numpy matrix
# MES: 2D list containing data from measurements.mes file
def buildP(MES, sigma0):
    n = len(MES)
    P = np.zeros([n, n])
    for i in range(0, n):
        P[i][i] = sigma0 / pow(MES[i][4], 2)

    return P


# Function to find coordinates in CNT by Point name
# Returns x, y as a tuple of floats
# CNT: 2D list containing data from the coordinates.cnt file
# name: name of point in Point column of CNT
def findPoint(CNT, name):
    P = None
    for point in CNT: # for all points
        if point[0] == name: # if point found
            P = Point(point)
    if not P: # if point could not be found
        exception_text = "Could not find point '" + name + "' in CNT"
        raise Exception(exception_text)
    return P


# Function to make the angles from atan2 between 0 and 2pi
# returns angle + npi in radians
# angle: input angle in radians
def npi(angle):
    while angle < 0:
        angle = angle + math.pi
    while angle > 2 * math.pi:
        estimated = angle - math.pi
    return angle

# Function to build the design and misclosure matrices
# returns A, w as a tuple of numpy matrices
# CNT: 2D list containing data from the coordinates.cnt file
# MES: 2D list containing data from measurements.mes file
# x: numpy array containing initial values of unknowns
# P: numpy weight matrix
def buildAw(CNT, MES, x):
    n = len(MES)  # n = number of measurements
    u = len(x)  # u = number of unknowns
    # initialize A and w to the proper size's
    A = np.zeros([n, u])
    w = np.zeros([n, 1])

    # loop through all measurements (all rows of A and w)
    for i in range(0,n):
        # get measurement variables
        mesID, mesInfo, mesType, mesValue = MES[i][0:4]

        # if angle measurement
        if mesType == 'Angle':
            # parse mesInfo
            try:
                Pto, Pat, Pfrom = mesInfo.split('_')
            except:
                exception_text = "Could not parse measurement info for ID = " + mesID
                raise Exception(exception_text)
            # get point coordinates
            Pj = findPoint(CNT, Pto)
            Pi = findPoint(CNT, Pat)
            Pk = findPoint(CNT, Pfrom)

            # for each point which in an unknown, replace (x,y) coordinates with updated values from unknowns vector x
            for P in [Pj, Pi, Pk]:
                if P.isUnknown():
                    x_index = 2*P.unknownNum(CNT)
                    P.x = x[x_index]
                    P.y = x[x_index + 1]

            # calculate row of A~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # partial derivatives:
            # if Pi is an unknown
            if Pi.isUnknown():
                # find which columns of the A matrix it should be in using the order of unknowns in CNT
                col = 2*(Pi.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pi
                A[i][col] = (Pk.y - Pi.y)/(pow(Pk.x-Pi.x, 2) + pow(Pk.y-Pi.y, 2)) - (Pj.y - Pi.y)/(pow(Pj.x-Pi.x, 2) + pow(Pj.y-Pi.y, 2))  # f/dxi
                A[i][col+1] = -(Pk.x - Pi.x)/(pow(Pk.x-Pi.x, 2) + pow(Pk.y-Pi.y, 2)) + (Pj.x - Pi.x)/(pow(Pj.x-Pi.x, 2) + pow(Pj.y-Pi.y, 2))  # f/dyi
            # if Pj is an unknown
            if Pj.isUnknown():
                col = 2*(Pj.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pj
                A[i][col] = (Pj.y - Pi.y)/(pow(Pj.x - Pi.x, 2) + pow(Pj.y - Pi.y, 2))  # f/dxj
                A[i][col+1] = -(Pj.x - Pi.x)/(pow(Pj.x - Pi.x, 2) + pow(Pj.y - Pi.y, 2))  # f/dyj
            # if Pk is an unknown
            if Pk.isUnknown():
                col = 2 * (Pk.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pk
                A[i][col] = -(Pk.y - Pi.y) / (pow(Pk.x - Pi.x, 2) + pow(Pk.y - Pi.y, 2))  # f/dxk
                A[i][col + 1] = (Pk.x - Pi.x) / (pow(Pk.x - Pi.x, 2) + pow(Pk.y - Pi.y, 2))  # f/dyk

            # calculate row of w~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # w = estimated - measured
            estimated = (math.atan2(Pk.y - Pi.y, Pk.x - Pi.x) - math.atan2(Pj.y - Pi.y, Pj.x - Pi.x))  # estimated angle in radians
            # make sure estimated angle is between 0 and 2pi
            estimated = npi(estimated)
            w[i] = estimated - mesValue

        # if distance measurement
        elif mesType == 'Dist':
            # parse mesInfo
            try:
                Pstart, Pend = mesInfo.split('_')
            except:
                exception_text = "Could not parse measurement info for ID = " + mesID
                raise Exception(exception_text)
            # get point coordinates
            Pi = findPoint(CNT, Pstart)
            Pj = findPoint(CNT, Pend)

            # for each point which in an unknown, replace (x,y) coordinates with updated values from unknowns vector x
            for P in [Pj, Pi]:
                if P.isUnknown():
                    x_index = 2*P.unknownNum(CNT)
                    P.x = x[x_index]
                    P.y = x[x_index + 1]

            # calculate row of A~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # partial derivatives:
            # if Pi is an unknown
            if Pi.isUnknown():
                # find which columns of the A matrix it should be in using the order of unknowns in CNT
                col = 2 * (Pi.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pi
                A[i][col] = (Pi.x - Pj.x)/pow(pow(Pi.x - Pj.x, 2) + pow(Pi.y - Pj.y, 2), 0.5)  # f/dxi
                A[i][col + 1] = (Pi.y - Pj.y) / pow(pow(Pi.x - Pj.x, 2) + pow(Pi.y - Pj.y, 2), 0.5)  # f/dyi
            # if Pj is an unknown
            if Pj.isUnknown():
                col = 2 * (Pj.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pi
                A[i][col] = (Pj.x - Pi.x) / pow(pow(Pi.x - Pj.x, 2) + pow(Pi.y - Pj.y, 2), 0.5)  # f/dxj
                A[i][col + 1] = (Pj.y - Pi.y) / pow(pow(Pi.x - Pj.x, 2) + pow(Pi.y - Pj.y, 2), 0.5)  # f/dyj

            # calculate row of w~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            estimated = pow(pow(Pi.x - Pj.x, 2) + pow(Pi.y - Pj.y, 2), 0.5)
            w[i] = estimated - mesValue
        else:
            exception_text = "Invalid measurement type for ID = " + mesID
            raise Exception(exception_text)

    return A, w

# Function to count the number of digits before the decimal place
# returns counts as an int if only 1 number is given, or as a list for multiple numbers
# nums: either a single int/float, or a list of ints/floats
def numlen(nums):
    try:  # if nums has no length (it's a float or int)
        tmp = len(nums)
    except:  # put it in an array
        nums = [nums]

    counts = []
    for num in nums:
        count = 0
        num = abs(num)
        while num >= 1:
            count = count + 1
            num = num/10.0
        counts.append(count)
    # if only 1 number was counted
    if len(counts) == 1:
        return counts[0]  # return as int
    else:
        return counts  # return as list


# Function to plot all coordinates from CNT
# returns fig and ax, the matplotlib figure and axis objects
# CNT: coordinates 2D list
# x: array with estimated unknown values
# fignumber: figure number that should be used
def plotCNT(CNT, x):
    unknownCol = '#fc4c4c'
    knownCol = '#03c2fc'
    ax = plt.axes()
    ax.axis('square')
    # save min/max x and y for plot limits
    # initialize to first point's (x,y)
    for i in range(0, len(CNT)):
        P = Point(CNT[i])
        c = knownCol  # set color to blue
        # if P is unknown
        if P.isUnknown():
            # set color to red if P is an unknown point
            c = unknownCol
            # update P(x,y) from unknowns vector x
            x_index = 2 * P.unknownNum(CNT)
            P.x = x[x_index]
            P.y = x[x_index + 1]
        # plot point
        plt.scatter(P.x, P.y, c=c)
        # add text
        plt.text(P.x+30, P.y+30, P.name)
    # set title
    ax.set_title('Adjusted Geodetic Network with Exaggerated Error Ellipses')
    # get min/max x
    minx = min([i[2] for i in CNT])
    maxx = max([i[2] for i in CNT])
    # get min/max y
    miny = min([i[3] for i in CNT])
    maxy = max([i[3] for i in CNT])
    # set limits with 10% buffer
    bufx = abs(maxx - minx)*0.1
    bufy = abs(maxy - miny)*0.1
    plt.xlim(minx - bufx, maxx + bufx)
    plt.ylim(miny - bufy, maxy + bufy)
    #plt.show(block=False)  # block=False allows the rest of the program to run while the plot is open
    return ax


# Function to draw an error ellipse on a figure
# returns nothing
# fig: matplotlib.pyplot.figure object
# x and y: (x,y) center of the ellipse
# semi_major: size of the semi-major axis
# semi_minor: size of the semi-minor axis
# theta: orientation of the semi-major axis in radians (angle counter-clockwise from the x axis to the semi-major axis)
# scale (optional): scale for the semi_major and semi_minor values to make the ellipse larger or smaller
# name (optional): set name for unknown point (used for zoomed in plot)
def drawEE(ax, x, y, semi_major, semi_minor, theta, scale=1, name=''):
    # convert theta to degrees
    thetad = theta*180/math.pi
    # create exaggerated ellipse on main plot:
        # ellipse is centered on the given (x,y)
        # width is the 2*semi-major axis
        # height in 2*semi-minor
        # angle takes degrees and rotates the ellipse counter-clockwise
        # theta is the counter-clockwise angle from the x axis, so we pass theta in degrees
    e = pat.Ellipse(xy=(x, y), width=2*semi_major*scale, height=2*semi_minor*scale,
                    angle=thetad, facecolor='none', edgecolor='red', label='Error Ellipse (exaggerated)')
    ax.add_artist(e)
    # create semi-major axis line
    l = Line2D([x, x + scale*semi_major*math.cos(theta)], [y, y + scale*semi_major*math.sin(theta)], color='red')
    ax.add_artist(l)

    # create zoomed in plot of unknown point with error ellipse to scale
    fig = plt.figure()
    ax2 = plt.axes()
    ax2.axis('square')
    plt.scatter(x, y)
    fig.suptitle("To-scale Error Ellipse for " + name + " (square plot)")
    e = pat.Ellipse(xy=(x, y), width=2*semi_major, height=2*semi_minor,
                    angle=thetad, facecolor='none', edgecolor='red')
    ax2.add_artist(e)
    # create semi-major axis line
    l = Line2D([x, x + semi_major*math.cos(theta)], [y, y + semi_major*math.sin(theta)], color='red')
    ax2.add_artist(l)
    # adjust limits to 2.2x the semi-major axis
    plt.xlim(x-semi_major*1.1, x+semi_major*1.1)
    plt.ylim(y-semi_major*1.1, y+semi_major*1.1)

    return


# Function to save all open figures as PDF
# I got this function from stack overflow user farenorth: https://stackoverflow.com/questions/26368876/saving-all-open-matplotlib-figures-in-one-file-at-once
# filename: filename of the pdf to be saved
# figs (optional): by default the function will grab all open figures. This option lets the user specify figures instead
def SaveFigs(filename, figs=None):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


# Function to find a single file of a certain filetype
# If more than 1 file is found, error is given
# returns a single filename
# ext: extension to look for (without '.')
def FindFiles(ext):
    files = os.listdir()
    found = []
    for file in files:
        if ext == file.split('.')[-1]:
            found.append(file)
    if len(found) > 1:
        raise Exception("More than 1 ." + ext + " file in directory. Either remove all but 1 ." + ext + """ file, or specify which file to use:
main(CNTFile='filename.cnt', MESFile='filename.mes')""")
    if len(found) == 0:
        raise Exception("No ." + ext + " file in directory. Ensure that ." + ext + """ file is in the root directory, or specify a path:
main(CNTFile='my/path/filename.cnt', MESFile='my/path/filename.mes')""")

    return found[0]  # return 1st string in found


#Function to test is the results are acceptable:
# rhat: residuals vector
# Crhat: variance covariance matrix for the residuals
# return results as a boolean list (true/false for each observation)
def sTest(rhat, Crhat):
    n = len(rhat)
    # check shape of Crhat
    if Crhat.shape[0] != n or Crhat.shape[1] != n:
        raise Exception("rhat must be nx1, and Crhat must be nxn")
    results = []
    # do test on each observation residual
    for i in range(0, n):
        if rhat[i] > 3*math.sqrt(Crhat[i][i]):
            results.append(False)
        else:
            results.append(True)
    if len(results) != n:
        raise Exception("results is the wrong size")
    # return results as a boolean list
    return results


# Function to write results from STest to file
# results: boolean list
# MES: measurements array (used to get measurement IDs)
# file: file object to output to
# returns nothing
def writeResults(results, MES, file):
    alltrue = True
    count = 0
    # loop through results
    for i in results:
        # if an observation failed
        if i == False:
            # change flag
            alltrue = False
            # get ID from MES
            ID = MES[count][0]
            # output error to file
            file.write("Measurement ID " + str(ID) + " failed\n")
        count = count + 1
    # if all observations passed
    if alltrue:
        # just write 1 line that tells the user all observations passed
        file.write("All observations passed!")
    return

