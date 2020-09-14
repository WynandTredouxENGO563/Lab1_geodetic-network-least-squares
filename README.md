# ENGO563 Lab 1
**Objective**: Conduct least squares adjustment on a geodetic network and analyze the results
## .cnt file
Contains all coordinates for points in the geodetic network. 

Format: `Point	Type	X	Y`

* **Point** is the point name or ID
* **Type** indicates if the point is a control point (C) or an unknown/estimated point (U)
* **X** and **Y** are the X and Y coordinates in meters

## .mes file
Contains all measurments

Format: `ID	Info	Type	Value	Std`
* **ID** is the measurement ID
* **Info** indicates where in the network the measurement was taken (A_P1_B indicates an angle measured at P1 from A to B, P1_A indicates a distance measurement from P1 to A)
* **Type** indicates what type of measurement it is (angle or distance)
* **Value** is the measurement value (DMS for angles seperated by a space, meters for distances)
* **Std** is the measurement standard deviation (decimal degrees for angles, meters for distances)
