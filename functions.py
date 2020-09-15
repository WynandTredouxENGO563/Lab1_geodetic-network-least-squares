# Function to read in text files for this project
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
                    if elements[2]=="Angle":
                        DMS = elements[3].split(' ')
                        elements[3] = float(DMS[0]) + float(DMS[1])/60 + float(DMS[2])/3600

                    row = [elements[0], elements[1], elements[2], elements[3], float(elements[4])]
                except:
                    raise Exception(exception_text)
                # save row in output
                output.append(row)

    return output