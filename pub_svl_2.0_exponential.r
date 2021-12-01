# This tool presents a new method to identify the vent location of tephra fall deposits 
# based on thickness or maximum clast size measurements.
# It is temporarily named "svl", which is short for Source Vent Locator.
# This file includes the associated R functions and their notations.

###
# Since the presented method gives estimate on the vent location and prevailing wind 
# direction, it can be used to produce the input for our method to model the thickness
# or maximum clast size distribution of 
# tephra fall deposits (Yang and Bursik, 2016; https://vhub.org/resources/3957).
###
# Author: Qingyuan Yang, Marcus Bursik, and E. Bruce Pitman.
###
# License: GPL. Use at your own risk. 
###
# The work was supported by National Science Foundation Hazard SEES grant number 1521855 to
# G. Valentine, M. Bursik, E.B. Pitman and A.K. Patra, and National Science Foundation
# DMS grant number 1621853 to A.K. Patra, M. Bursik, and E.B. Pitman.
###
# We appreciate your comments, suggestions, and feedback. 
# Please feel free to contact us through vhub or email: qyang5@buffalo.edu  
# Notes:
#	The functions below are the version of the method that is coupled with the semi-	
#	empirical model proposed by Yang and Bursik (2016).

###### Functions ######
#################################Intermediate functions#################################
#################### Function name: sp_par_cal ####################
# Function to calculate downwind and cross wind distance of the input sample sites given known
# source vent and wind direction.
# The name "sp_par_cal" is short for "spatial parameter calculation".
# Input:
	# obs: 
	#	n-by-3 matrix of x, y, z, indicating coordinates and measured thickness 
	#	or maximum clast size of each sample site. x and y should be in utm coordinates,
	#	and z should be in millimeter with a minimum value of 1 mm.
	# sv: 
	#	a vector of two elements (utm coordinates) indicating the assumed vent location.
	# wind_dir
	# 	a numeric value of assumed wind direction (from north clockwise).
# Output:
	# an n-by-7 matrix, 
	# The 1st and 2nd columns: 
	#	original coordinates of sample sites;
	# The 3-5th columns:
	#	total distance (dist), downwind (dd), and crosswind (cd) distance 
	# 	with respect to the source vent (sv).
	# The 6th column:
	#	difference between total distance and downwind distance.  
	# The 7th column:
	#	thickness or grain size measurement at each sample site.
	###### Brief description ######
	# This function uses basic vector/matrix multiplication in R to prepare the data for
	# calculating and fitting the objective function, GIVEN ASSUMED FIXED source 
	# vent location and wind direction.
sp_par_cal <- function(obs, sv, wind_dir){
	sv = as.vector(sv)
	zr <- obs[,3]
	xy <- as.matrix(obs[,c(1,2)])
	nr <- nrow(obs)
	sv_matrix <- matrix(ncol=2, nrow=nr)
	sv_matrix[,1] <- rep(sv[1], nr)
	sv_matrix[,2] <- rep(sv[2], nr)
	unit_vector <- as.matrix(c(sin(wind_dir*pi/180),cos(wind_dir*pi/180)))
	rel_xy <- xy - sv_matrix
	dist <- (rel_xy[,1]^2 + rel_xy[,2]^2)^0.5
	dd <- (rel_xy) %*% (unit_vector)
	cd <- (dist^2 - dd^2)^0.5     
	output <- matrix(nrow = nr, ncol=7)
	output[,1:2] <- xy
	output[,3:5] <- cbind(dist, dd, cd)
	output[,6] <- dist - dd  
	output[,7] <- zr
	return(output)
}

#################### Function name: fi ####################
# Function to calculate the sum of squared residuals (calculated under log-scale) 
# GIVEN FIXED source vent location and wind direction.	
# The name "fi" is short for "fitting".												
# Input:
	# obs, sv, and wind_dir: 
	#	inputs for function "sp_par_cal".
	###### Output ######
	# ssr:
	#	a numeric value which is the sum of squared residuals (ssr) from the fitting.
	# Note:
	#	The method uses the semi-empirical model of Yang and Bursik (2016).
	#	Note the line with "***" on the right, which is used for the
	#	semi-empirical model proposed by Gonzalez-Mellado and De la Cruz-Reyna 
	# 	(2010). It is not used here.
fi <- function(obs, sv, wind_dir){
	td <- as.data.frame(sp_par_cal(obs, sv, wind_dir))		
		# Call the "sp_par_cal" function, and turn the output to data.frame. 	
	colnames(td) <- c("x","y","dist","dd","cd", "dif","zr") 	
		# Name the columns.
	td$z = log10(td$zr)						
		# Transform the thickness into log-scale.
	## Fitting and calculate the ssr
#       fit.res <- lm(z~dif+log(dist),data=td)		# *** Case 1: power-law model, not used here.
        fit.res <- lm(z~dist+dd,data=td)		# Case 2: exponential model.
      	ssr <- sum(fit.res$residuals^2)			# Calculate the ssr.
	return(ssr)
}

#################### Function name: fi2 ####################							
# This function is basically identical to the function "fi". 
# The only difference is the output.
# It is called when the vent position is found. 
# In addition to "ssr" in function "fi", it also returns the fitted coefficients 
# of the semi-empirical model and the associated r-square value.
# With these values, we could have a better understanding on the predicted result.
# Input
	# Identical to the inputs for function "fi"
# Output
	# a vector with five elements including:
	# ssr: sum of squared residuals;
	# fit.res$coefficients: fitted coefficients (three values) of the semi-empirical model given
	# 	fixed vent location and wind direction;
	# summary(fit.res)$r.squared): r-squared value of the fitted result given
	#	fixed vent location and wind direction (a numeric value).
fi2 <- function(obs, sv, wind_dir){	
	td <- as.data.frame(sp_par_cal(obs, sv, wind_dir))
	colnames(td) <- c("x","y","dist","dd","cd", "dif","zr")
	td$z = log10(td$zr)
	## Fitting and measure the ssr
#      	fit.res <- lm(z~dif+log(dist),data=td)		# Case 1: power-law model, not used here.
       	fit.res <- lm(z~dist+dd,data=td)		# Case 2: exponential model.
      	ssr <- sum(fit.res$residuals^2)
	return(c(ssr, fit.res$coefficients, summary(fit.res)$r.squared))
}


#################### Function name: cg ####################
# Function to generate points (in the x-y plane) around a given point. This given point
# is the proposed vent location from the previous iteration.
#
# The generated points are proposed such that the following function compares if 
# one of them or the proposed vent location from the previous iteration is closer to the true vent location. 
# The generated points are located in the four cardinal directions with respect to the proposed vent location
# from the previous iteration. 
# The function name is short for "Circle Generator".  
# Input:
	# sv: 
	#	the proposed vent location from previous iteration.
	#	Its form is identical to "sv" in "sp_par_cal".
	# h:
	#	a numeric value (search radius), indicating how far the four generated points are 
	#	from the previous proposed vent location. 

# Output
	#	a 4-by-2 matrix containing the coordinates of the four generated points.
cg <- function(sv, h){					
	output  = cbind(c(rep(sv[,1],4)+sin(0:3*90*pi/180)*h),
		    c(rep(sv[,2],4)+cos(0:3*90*pi/180)*h))
	return(output)
}

#################### Function name: baf_pre ####################
# This function prepares for the function "baf" below. It is designed 
# to avoid local minima in estimating the wind direction.
# The name "baf_pre" is short for "best-fitted angle finder_prepare".
# Input
	# sv and obs:
	#	inputs required for function "sp_par_cal" 
	# nlps:
	#	A numeric value indicating the number of iterations
	# 	used to estimate the prevailing wind direction.
	# 	It is not used here, but is necessary for the ongoing functions.
# Output
	# A vector containing two or three values:
	# 	If the first element is 1, 
	# 		no local minima occur. The 2nd element is a rough estimate
	#		on the prevailing wind direction.
	#	If the first element is 2,
	# 		local minima are present. The next two elements are rough estimates
	#		on the global and local minima of wind direction. Both of them will be applied to ongoing
	#		functions.
baf_pre <- function(sv, obs, nlps){
	pot_ang <- as.data.frame((1:36)*10)									
	pot_ssr <- apply(pot_ang, MARGIN = 1, FUN=fi, obs = obs, sv = sv)
	#lowest_rank = which(rank(pot_ssr)==1)			# Power-law ***
	#next_lowest_rank = which(rank(pot_ssr)==2)		# Power-law ***
	lowest_rank = min(which(rank(pot_ssr)<2))				#exponential
	next_lowest_rank = min(which(rank(pot_ssr)<4 & rank(pot_ssr)>2))	#exponential
	if(abs(lowest_rank - next_lowest_rank)==1 || abs(lowest_rank - next_lowest_rank)==35){
		return(c(1, pot_ang[lowest_rank,1]))
	}else{
		return(c(2, pot_ang[lowest_rank,1], pot_ang[next_lowest_rank, 1]))
	}
}

#################### Function name: baf_processor ####################
# Function to estimate the prevailing wind direction with an assumed and
# fixed vent location. It uses a standard one-dimensional gradient descent method.
# The name "baf_processor" is short for "best-angle finder_processor".
# Input
	# sv and obs: inputs for the function "sp_par_cal".
	# nlps: number of loops for the 1d gradient descent method.
	# 	It is sufficient to set its value to 20.
	# proposed_winddir: the proposed wind direction from the function "baf_pre".
	# This function uses the output from function "baf_pre" as 
	# the initial guess on wind direction.
# Output
	# a 1-by-2 matrix containing: 
	#	the predicted wind direction;
	#	the corresponding sum of squared residuals.
baf_processor <- function(sv, obs, nlps, proposed_winddir){										
	fdm <- matrix(c(-1,0,1,-1,0,1), nrow=2)					
	h <- 5				  
	warns <- 0					
	resm <- data.frame(matrix(ncol=2,nrow=2))
	gdnt <- c()	
	wind_dir = proposed_winddir
	res <- fi(obs, sv, wind_dir)
	resm[1,1] <- wind_dir
	resm[1,2] <- res[1]
	wind_dir <- wind_dir + h
	res <- fi(obs, sv, wind_dir)
	resm[2,1] <- wind_dir
	resm[2,2] <- res[1]
	if(resm[1,2] >= resm[2,2]){
		lid <- 1	
	}else{
		lid <- -1
	}
	r = 1
	# lid=1 -> going forward, wind_direction+ ;
	# lid=-1 -> going backward, wind_direction- ;
	# lid=0 -> going in between, the wind direction is between the last two values;
	# r=0 -> h has just been shrunk;
	# r=1 -> the h has been shrunk for the last two iterations.
	for(i in 3:nlps){
		if(lid == 1){
			wind_dir <- resm[order(resm[,2]),][1,1] + h
			lid <- 1
		}else if(lid == -1){
			wind_dir <- resm[order(resm[,2]),][1,1] - 2 * h
			lid <- -1
		}else if(lid == 0 && r == 1){
			h <- h * 0.4
			wind_dir <-resm[order(resm[,2]),][1,1] + h
			r = 0
		}else if(lid == 0 && r == 0){
			wind_dir <-resm[order(resm[,2]),][1,1] - h
			r = 1
		}
	
		res <- fi(obs, sv, wind_dir)
		resm[i,1] <- wind_dir
		resm[i,2] <- res[1]

		resm <- resm[order(resm[,1]),]
			if(lid == 1){
				gdnt <-  fdm %*% resm[c(i-2,i-1,i),2]
			}
			if(lid == -1){
				gdnt <- fdm %*% resm[c(1,2,3),2] 
			}
			if(lid == 0 ){	
				n <- which(resm[,2]==min(resm[,2]))
				gdnt <- fdm %*% resm[c(n-1,n, n+1) ,2]
			}
		if(gdnt[1]<0 && gdnt[2]>=0){lid <- 0}
		if(gdnt[1]>=0 && gdnt[2] >= 0){lid <- -1}
		if(gdnt[1]<0 && gdnt[2] <0){lid <- 1}
		if(gdnt[1] > 0 && gdnt[2] < 0){lid <- -1}
		
		if(h < 0.01){break}
	}
	temp <- resm[which(resm[,2]==min(resm[,2])),]
	output <- as.matrix(cbind(temp[1], temp[2]), nrow = 1)
	
	return(output)		
}	


#################### Function name: baf ####################
# As noted above, local minimum could occur for estimating the prevailing wind direciton.
# Therefore this function is designed to receive outputs from the function "baf_pre".
# It first identifies if local minima occur or not, then applies the proposed initial
# guess(es) to the function "baf_processor", and collects the output.
# If local minima occur, it compares the resultant two sum of squared residuals to
# determine which one represents the true prevailing wind direction (global minimum).  
# Input
	# sv and obs:
	#	inputs required for function "sp_par_cal".
	# nlps:
	#	a numerical value required for the function "baf_processor",
	#	as described before. 
# Output
	# Optimum output from function "baf_processor".
baf <- function(sv, obs, nlps){
	raw_inference <- baf_pre(sv, obs, nlps)
	if(raw_inference[1] == 1){
		output = baf_processor(sv, obs, nlps, raw_inference[2])
		return(output)
	}else{
		output_1 = baf_processor(sv, obs, nlps, raw_inference[2])
		output_2 = baf_processor(sv, obs, nlps, raw_inference[3])
		if(output_1[2]<output_2[2]){return(output_1)}else{return(output_2)}
	}
}

#################### Function name: compare ####################
# This function COMPAREs if the proposed vent location from the last iteration
# or one of the points generated from function "cg" is closer to the true vent location.	
# 	If the proposed vent location from the last iteration is closer to the true vent,
# 	it passes this proposed vent location to the next iteration.
# 	If one or more points generated from function "cg" are closer
# 	to the true vent, it provides a direction (unit vector) which is used to propose 
#	the vent location for the next iteration (which will lie in the direction represented
#	by the unit vector).
# Input
	# sv and obs:
	#	inputs required for function "sp_par_cal".
	# h:
	# 	input required for function "cg".
	# nlps:
	#	inputs required for functions "baf_pre", "baf_processor", and "baf".
	# h_ratio:
	#	a numeric value that can be used to further constrain the search radius in function "cg".
	#	From our experiments, it seems that to set this value to 1 is an appropriate option. 
# Output
	# 	if: 
	#	the proposed vent location from the last iteration is closer to the true vent:
	#		the output is a numeric value 1. We keep this point for the next iteration.
	#	else:
	#		the output is a vector of three elements. The first element is the numeric value 2.
	#		The second and third elements together form a unit vector, and the proposed 
	#		vent location for the next iteration will be at that direction with respect to 
	#		the current proposed vent location.  
	#		The way the unit vector is determined can be thought of as a crude approximation 
	#		of the divergence of the cost function.
compare <- function(sv, obs, h, nlps, h_ratio){
	rsm <- (cg(sv, h*h_ratio))					
	cr <- baf(sv, obs, nlps)				

	rr = t(apply(X = rsm,  MARGIN = 1, FUN = baf, obs = obs, nlps = nlps))
	first_order_grad <- rr[,2]-cr[2]

	if(min(first_order_grad>0)){
		return(1)
	}else{
		first_order_grad[first_order_grad > 0] = 0	
		direction = c(1,1,-1,-1)
		grad_y = (-1*direction[c(1,3)])%*%first_order_grad[c(1,3)]
		grad_x = (-1*direction[c(2,4)])%*%first_order_grad[c(2,4)]
		grad_unit = c(grad_x/(grad_x^2+grad_y^2)^0.5, grad_y/(grad_x^2+grad_y^2)^0.5)
		return(c(2,grad_unit))
	}
}
####################End of intermediate functions#################################

#################### Function name: gd_simplified ####################
# This function integrates the above functions, implements ANOTHER gradient descent 
# method in the x-y plane, and estimates the source vent location.   
# See notations in the file "demo_svl_2.0.R" for the description of this 
# function, and how to set values for the input parameters.
gd_simplified <- function(sv, obs, runs, h, r, nlps, h_ratio, numb){		
 
	obs = obs[numb,] 
	tcoord <- as.data.frame(matrix(ncol=2))			
	dec <- c()
	hl <- c()
	ihd <- h

	ir=1							
	jg <- 0

	warns <- 0
	for(i in 1:runs){
			dec <- compare(sv, obs, h, nlps, h_ratio)
			if(dec[1]!= 1 && jg ==0){
					sv <- sv + h*dec[2:3]
					tcoord[i,] <- sv
					hl[i] <- h					
					jg <- 0
			}else if(dec[1] != 1 && jg == 1){
					h <- h *(1 - r)
					sv <- sv + h*dec[2:3]
					tcoord[i,] <- sv
					hl[i] <- h					
					jg <- 0
			}else if(dec[1] == 1 && h/ihd >= 0.0001){	
				h = h * r
				hl[i] <- h	
				jg <- 1
			}	
		if(h < 0.1){break}
	}
	bd <- baf(sv, obs, nlps)
	rownames(tcoord) <- NULL

	output_sv = tcoord[nrow(tcoord),]
	output_ang = bd[1]

	fitted_results= fi2(obs, as.matrix(sv), output_ang)

	return(as.matrix(cbind(sv, output_ang, h, coef = t(fitted_results[c(2,3,4)]), ssr = fitted_results[1], rsquare = fitted_results[5])))
}
########################################## End of functions #####################################################

# References:
#	Gonzalez-Mellado, A. O., and S. De la Cruz-Reyna. "A simple semi-empirical approach to model thickness of 
# ash-deposits for different eruption scenarios." Natural Hazards and Earth System Sciences 10.11 (2010): 2241.
# 	Yang Q, Bursik M. A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall 
# deposits. Bulletin of Volcanology, 2016, 78(10): 68.

