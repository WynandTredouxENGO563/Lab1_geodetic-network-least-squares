import numpy as np


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
        for line in lines[header_lines:-1]:
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

                    row = [elements[0], elements[1], elements[2], elements[3], float(elements[4])]
                except:
                    raise Exception(exception_text)
                # save row in output
                output.append(row)

    return output


# Function to build the unknowns vector x
# returns x as a numpy array
# CNT: 2D list containing data from the coordinates.cnt files
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


# Function to build the design and misclosure matrices
# returns A, w as a tuple of numpy matrices
# CNT: 2D list containing data from the coordinates.cnt files
# MES: 2D list containing data from measurements.mes file
# x: numpy array containing initial values of unknowns
# P: numpy weight matrix
def buildAw(CNT, MES, x, P):
    n = len(MES)  # n = number of measurements
    u = len(x)  # u = number of unknowns
    # initialize A and w to the proper size's
    A = np.zeros([n, u])
    w = np.zeros([n, 1])
    # loop through all measurements
    for i in range(0,n):
        # get measurement
        mesType = MES[i,]

    return A, w
