from classes import *
import numpy as np
import math


# Function to read in text files for this project
# returns either a string (for .txt files) or a 2D list (for .cnt or .mes files)
# Filesname: name of file including extension
# header_lines: number of lines that the header takes up (will be skipped)
def readfile(filename, header_lines):
    # acceptable file extensions
    ext = ['txt', 'mes', 'cnt']
    # get file extension
    tmp = filename.split('.');
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
                    # convert DMS to DD for angle measurements
                    if elements[2] == "Angle":
                        DMS = elements[3].split(' ')
                        elements[3] = float(DMS[0]) + float(DMS[1]) / 60 + float(DMS[2]) / 3600

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
            P = Point(point[0], point[1], point[2], point[3])
    if not P: # if point could not be found
        exception_text = "Could not find point '" + name + "' in CNT"
        raise Exception(exception_text)
    return P

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
        mesID, mesInfo, mesType, mesValue, mesStd = MES[i][0:5]

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
                A[i][col+1] = (Pk.x - Pi.x)/(pow(Pk.x-Pi.x, 2) + pow(Pk.y-Pi.y, 2)) - (Pj.x - Pi.x)/(pow(Pj.x-Pi.x, 2) + pow(Pj.y-Pi.y, 2))  # f/dyi
            # if Pj is an unknown
            if Pj.isUnknown():
                col = 2*(Pj.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pj
                A[i][col] = (Pj.y - Pi.y)/(pow(Pj.x - Pi.x, 2) + pow(Pj.y - Pi.y, 2))  # f/dxj
                A[i][col+1] = (Pj.x - Pi.x)/(pow(Pj.x - Pi.x, 2) + pow(Pj.y - Pi.y, 2))  # f/dyj
            # if Pk is an unknown
            if Pk.isUnknown():
                col = 2 * (Pk.unknownNum(CNT))
                # calculate partial derivative for (x,y) of Pk
                A[i][col] = (Pk.y - Pi.y) / (pow(Pk.x - Pi.x, 2) + pow(Pk.y - Pi.y, 2))  # f/dxk
                A[i][col + 1] = (Pk.x - Pi.x) / (pow(Pk.x - Pi.x, 2) + pow(Pk.y - Pi.y, 2))  # f/dyk

            # calculate row of w~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # w = estimated - measured
            estimated = (math.atan2(Pk.y - Pi.y, Pk.x - Pi.x) - math.atan2(Pj.y - Pi.y, Pj.x - Pi.x))*180/math.pi  # estimated angle in degrees
            # make sure estimated angle is between 0 and 360 degrees
            while estimated < 0:
                estimated = estimated + 180
            while estimated > 360:
                estimated = estimated - 180
            w[i] = estimated - MES[i][3]

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
            w[i] = estimated - MES[i][3]
        else:
            exception_text = "Invalid measurement type for ID = " + mesID
            raise Exception(exception_text)

    return A, w
