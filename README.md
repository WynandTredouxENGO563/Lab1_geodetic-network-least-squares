# ENGO563 Lab 1
**Objective**: Conduct least squares adjustment on a geodetic network and analyze the results

## Input files
### .cnt file
Contains all coordinates for points in the geodetic network. 

Format: `Point	Type	X	Y`

* **Point** is the point name or ID
* **Type** indicates if the point is a control point (C) or an unknown/estimated point (U)
* **X** and **Y** are the X and Y coordinates in meters
### .mes file
Contains all measurements

Format: `ID	Info	Type	Value	Std`
* **ID** is the measurement ID
* **Info** indicates where in the network the measurement was taken (A_P1_B indicates an angle measured at P1 from A to B, P1_A indicates a distance measurement from P1 to A)
* **Type** indicates what type of measurement it is (angle or distance)
* **Value** is the measurement value (DMS for angles separated by a space, meters for distances)
* **Std** is the measurement standard deviation (arc seconds for angles, meters for distances)

## Output Files
### output.out
This files gives all relevant info about the adjustment
It consists of 9 sections
* **Section 1** - Unnamed
	* Execution date and time
	* number of iterations
	* threshold value used (calculated as 1/2 the smallest measurement standard deviation)
	* sigma0 - a priori variance factor (1 by default)
* **Section 2** - Observations/Unknowns Summery
	* Number of measurements
	* Number of unknowns
	* Total degrees of freedom
	* Posteriori/unit Variance Factor
* **Section 3** -  Estimated Unknowns
	* Estimated coordinates of unknowns in meters
	* difference in x and y from the initial estimates given in .cnt (also in meters)
* **Section 4** - Vector of Residuals
	* Measurement residuals in radians and meters (corrected - original measurements)
* **Section 5** - Vector of Corrected Measurements
	* Corrected measurements in radians and meters
* **Section 6** - Variance Covariance Matrix of Unknowns
	* Variance/Covariance matrix of unknowns in m<sup>2</sup> 
* **Section 7** - Variance Covariance Matrix of Corrected Measurements
	* Variance/Covariance matrix of unknowns in radians<sup>2</sup> and m<sup>2</sup> along the diagonal
	* radians<sup>2</sup>, m<sup>2</sup>, or m*rads in the off-diagonal terms
* **Section 8** - Variance Covariance Matrix of Residuals
	* Same units as section 7
* **Section 9** - Error Ellipse
	* Semi-major and semi-minor axis lengths in meters
	* semi-major axis azimuth in radians (clockwise angle from N/y-axis to the semi-major axis)
