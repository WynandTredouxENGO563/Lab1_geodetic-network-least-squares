# Point class
# stores all info about a point from a row of CNT
class Point:
    # constructor: pass list of [name, type, x, y]
    def __init__(self, list):
        try:
            self.name = list[0]
            self.type = list[1]
            self.x = list[2]
            self.y = list[3]
        except:
            raise Exception()

    # Function to check if this point in an unknown
    def isUnknown(self):
        if self.type == 'U' or self.type == 'u':
            return True
        else:
            return False

    # Function to check which unknown this is
    def unknownNum(self, CNT):
        # if self is not an unknown
        if not self.isUnknown():
            # give error
            exception_text = "Point '" + self.name + "' is not an unknown"
            raise Exception(exception_text)
        # loop through CNT and count unknowns until self is reached
        count = 0
        for point in CNT:  # for all points in CNT
            P = Point(point)
            if P.name == self.name:  # if point is reached in CNT
                return count # return the count of unknowns
            elif P.isUnknown(): # if point in CNT is an unknown
                count = count + 1 # increment

        # if all points in CNT were looped through and self was not found
        # error
        exception_text = "Point '" + self.name + "' can not be found in CNT"
        raise Exception(exception_text)