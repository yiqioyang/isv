# This file includes a demonstration of how to use the method "svl" to identify the source
# vent location of tephra fall deposits based on thickness or maximum clast size
# measurements. 
# Author: Qingyuan Yang, Marcus Bursik, and E. Bruce Pitman.
# GPL: License. Use at your own risk. 
###
# The work was supported by National Science Foundation Hazard SEES grant number 1521855 
# to G. Valentine, M. Bursik, E.B. Pitman and A.K. Patra, and National Science Foundation
# DMS grant number 1621853 to A.K. Patra, M. Bursik, and E.B. Pitman.
# We appreciate your comments, suggestions, and feedback. 
# Please feel free to contact us through vhub or email: qyang5@buffalo.edu  
########################################## Demonstration ################################
# In this demonstration, we use a thickness dataset of the North Mono Bed 2 
# (Sieh and Bursik, 1986) as an example. North Mono Bed 2 was erupted from Upper Dome, which is a
# part of the Mono-Inyo Craters (eastern central California).
# It is noted that the method can also be applied to maximum clast size measurements.

# Set working directory.
setwd("/the/directory/where/you/put/your/data/") 

#---

# The input dataset should have three columns:
#	The first two columns are the coordinates (x or E and y or N) of sample sites in the UTM system.
#	The third column should include the corresponding thickness or maximum clast size 
#	measurements in millimeters, with a minimum value of 1 mm.
# The input dataset should be a ".csv" file, and separated by ",".
# The head or column names of the three columns should be "x", "y", and "zr".
# In the case of  North Mono Bed 2, the x and y coordinates are in UTM Zone 11. 
# Read dataset from the working directory.      
td = read.table("nmb2.csv", header=TRUE, sep = ",")

#---      

# Check data structure (optional command).
str(td)
# Return:
	#'data.frame':      114 obs. of  3 variables:
	#$ x : num  320364 328807 327410 329283 327802 ...
	#$ y : num  4194955 4192135 4193500 4194283 4194929 ...
	#$ zr: num  0 0 0 0 0 0 0 0 0 0 ...
# NOTE: if the data include observations of zero thickness; it is necessary to exclude them for the 
#	method to work (optional command).
#
#
#---

# Exclude zero-thickness observations.
td = subset(td, zr!=0)
      
#---

# Check data structure in a different way (optional command).
# Command below shows the first six rows of the dataset.
head(td)
# Return:
	#          x       y zr
	#12 326341.2 4199522 22
	#15 327685.3 4203216 10
	#26 320851.0 4197091 38
	#27 321276.5 4197466 82
	#28 320946.3 4197606 96
	#29 321606.7 4197460 78
# Note: the 1st column with # is just sample numbers; these are not
# part of the dataset. The data are not in order (1, 2, 3,...)
# because zero-thickness observations have been excluded. 

#---

# Define the initial guess of source vent location.
# These coordinates should be within the area where the vent is likely 
# to be.  
# We recommend users to try different coordinates within that area to 
# run the method,
# and then check if the results converge to the same point.
sv_assumed = cbind(x = 322227, y = 4181034)
# In the example case, assuming that we known the vent is from the Mono-Inyo Craters,
# the values of the initial guess are the coordinates of Obsidian Dome, which is > 14 km
# south of Upper Dome, the true vent location of North Mono Bed 2.

#---

# Source functions of the method "svl".
# Note: the semi-empirical model proposed by Yang and Bursik (2016) 
# is used here.
source("/where/you/put/your/source/code/pub_svl_2.0_exponential.r")

#---

# Run the method "svl".

# Brief description of the method:
#	Starting with an initial guess, "sv_assumed", on vent location, the
#	value is updated in each iteration (loop).					  
# Start loop:								 <----------------|
#	The method proposes different possible vent locations around "sv_assumed"   	  |	
#	that are at a certain distance (distance: h = 1000 m in this case) from it,       |
#	and compares if they are closer or farther to the true vent compared 		  |
#	with "sv_assumed".								  |
#          										  |
#	If "sv_assumed" is closer to the true vent:                                       |
#		we shrink the search radius (controlled by r * h, 0.7 * 1000 = 700 here)  |
#		to improve the resolution,						  |
#		and do the same thing with the smaller search radius; --------------->----/
#	If one of the nearby points is closer to the true vent:				  |	
#		the guess on vent location, "sv_assumed", 				  |
#		is updated/replaced (see function "compare" for more details),		  |
#		and with the same search radius, the method does			  |
#		the same thing with the new guess on vent location.------------------>----/
#	
#	In the actual implementation of the method, the search radius is updated
#	in a slightly more efficient manner (see detailed description in our work, whose 
#	link will be provided soon). 
#	This does not affect the notation in this file.
#
#	As described above, after each iteration, 
#	the search radius gets either smaller or unchanged. 
#	When the search radius is sufficiently small, the vent is found.
# End loop.

#	The number of iterations is controlled by "runs".
#---
#       How to assign values to the input parameters:
#	sv = sv_assumed, the initial guess on the source vent location;
#	td: the input dataset;
#	runs: the number of iterations; 
#	h: the initial search radius (in meters);
#	r: the shrink ratio in the search radius. 
#		Suggested value: 0.7
#	nlps, h_ratio: these are related to estimating the wind direction in each
#		iteration. It is sufficient to keep them fixed as 20 and 1, respectively. 
#	numb: the number of rows of data points within "td" that will be used as input. 
#	      To use the complete dataset, set it to "1:nrow(td)".
#---
#***	Values of "runs" and "h" need to be selected with care. ***
#	Here we use the example of the North Mono Bed 2 to show how the values are determined.
#
# 	Assuming that we know the tephra was erupted from the Mono-Inyo Craters, 
#	which cover a range of >15 km in the y direction, and about 5 km in the x-direction.
#
#	By specifying the initial search radius to be 1000 m (h = 1000), we ensure that it takes at most 
#	~15 steps/iterations before the method shrinks the search radius. 
#	By the time it shrinks the search radius, it can be assumed approximately that the true vent
#	is within ~1500 m of the current point. 
#
# 	Consider the worst-case scenario, namely the search radius keeps shrinking.  This means that
#	the searching resolution is improving in each iteration at the lowest rate.
#	In this case, the initial search radius keeps shrinking at a rate of 0.7 (r = 0.7). 
#	For the resolution to be at ~10 m scale, the method requires 13 more iterations: 1000*(0.7)^13 = 9.61.
# 	Based on this, it is ensured that runs = 40 (greater than 15+13) 
#	iterations are sufficient.
#
# 	Based on the argument above, we recommend users to follow this strategy to assign values of
#	"h" and "runs":
#	h = (the maximum height or width of the region of interest)/10
#	For "runs", if we want to reach the solution at a ~10 m scale, then 
#	we need to find out the value "x" in the brackets of the function: h*0.7^(x) = ~10. 
#	And then set "runs = 10 + x + E", where E could be an integer ranging from 5 to 20, a fudge-factor. 
#	Greater value in E means more additional iterations, 
#	which is likely to yield a more accurate estimate. 
result_linear = gd_simplified(sv = sv_assumed, td, runs = 40, h = 1000, r = 0.7, nlps = 20, h_ratio = 1, numb = 1:nrow(td))

# Check the output.
result_linear
# Return:
      #       x       y output_ang          h (Intercept)          dist            dd      ssr   rsquare
      #322224.4 4197601   172.3507 0.06281024    2.231973 -0.0003055393 -0.0002597007 2.877154 0.9419849
      
# Interpretation:
#	x and y: 	estimated vent coordinates of the North Mono Bed 2. This location is within the extent of 
# 			its true source, Upper Dome.
#	output_ang: 	estimated wind direction (from north clockwise). In this case, it is the UPWIND direction.
#			See note below on how to tell if the estimated wind direction
#			is downwind or upwind.
#	h:		the length of the search radius (m) from the last iteration. If this value is small
#			enough (e.g., <10), then the method converged. 
#			If this value is comparable to the initial search radius (say 490 = 1000*0.7*0.7), 
#			then users should try to increase the number of iterations, 
#			namely increase the value "runs" (e.g., set it to be 60),
#			and run the method again.
#				If the method is run with more iterations, 
#				but the resultant "h" is still comparable 
#				to the initial search radius, that means the method fails
#				to converge. 
#				Then users should try to change the initial
#				guess on vent location (sv_assumed), and run the method.
#	(intercept):	fitted coefficients for the semi-empirical model.
#	dist:		same as above.	
#	dd:		same as above.
#	ssr:		sum of squared residuals from the fitting given estimated vent location and  
#			wind direction.
#	rsquare:	r-squared value for the final estimate given estimated vent location and wind direction.
# 
# Given sparse data, it is important to check if the fitted coefficients are physical.
# For the exponential model, this can be done by checking if dist+dd<0. If their sum dist+dd > 0 or dist > 0, this means that 
# along the dispersal direction, the thickness or maximum clast size increases with distance, which is generally unphysical.  

#------------------------------------------
# Use the semi-empirical (power-law) model proposed by Gonzalez-Mellado and De la Cruz-Reyna (2010).
# Source the functions of the method "svl".
source("/where/you/put/your/source/code/pub_svl_2.0_power-law.r")

#---

# Run the method "svl" with identical inputs.
result_power = gd_simplified(sv = sv_assumed, td, runs = 40, h = 1000, r = 0.7, nlps = 20, h_ratio = 1, numb = 1:nrow(td))
# Check the result.
result_power
# Return:
# 	x              y output_ang         h (Intercept)           dif  log(dist)      ssr   rsquare
# 	322859.7 4198078   337.9795 0.5338782    5.194675 -0.0001700573 -0.4008759 2.084604 0.9579659

# Interpretation:
# 	See lines 166-190.
#	The estimated vent location is within the area of Upper Dome, the true vent of 
#	this tephra deposit.
#	### the output_ang in this case is the DOWNwind direction. 
#	See note below on how to tell if the estimated wind direction is downwind or upwind.
# 	If log(dist) > 0 or (Intercept) < 0, then unphysical prediction occurs.

######

#---

# Note on how to tell if the estimated wind direction (output_ang) is downwind or upwind:
# Due to the design of "svl", the output "output_ang" could be either the upwind or downwind direction.
# This is related to the gradient descent method adopted in "svl". 
# To tell if it is upwind or downwind,
# users need to examine the relationship between the fitted coefficients:

	# If "svl" is combined with the exponential method of Yang and Bursik (2016):
	#	if dd < 0, then "output_ang" is the upwind direction;
	#	if dd > 0, then "output_ang" is the downwind direction.

	# If "svl" is combined with the power-law method of Gonzalez-Mellado and De la Cruz-Reyna (2010):
	#	if dif < 0, then "output_ang" is the downwind direction;
	#	if dif > 0, then "output_ang" is the upwind direction.

# References:
#	Gonzalez-Mellado, A. O., and S. De la Cruz-Reyna. "A simple semi-empirical approach to model thickness of 
# ash-deposits for different #eruption scenarios." Natural Hazards and Earth System Sciences 10.11 (2010): 2241.
# 	Sieh K, Bursik M. Most recent eruption of the Mono Craters, eastern central California. Journal of Geophysical 
# Research: Solid Earth, 1986, 91(B12): 12539-12571.
# 	Yang Q, Bursik M. A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall 
# deposits. Bulletin of Volcanology, 2016, 78(10): 68.
