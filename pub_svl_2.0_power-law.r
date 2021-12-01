# This offline tool is a new method to identify the vent location of tephra fall deposits 
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
#	empirical model proposed by Gonzalez-Mellado and De la Cruz-Reyna (2010).
###### Functions ######
#################################Intermediate functions#################################
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

fi <- function(obs, sv, wind_dir){														
	td <- as.data.frame(sp_par_cal(obs, sv, wind_dir))
	colnames(td) <- c("x","y","dist","dd","cd", "dif","zr")
	td$z = log10(td$zr)

      	fit.res <- lm(z~dif+log(dist),data=td)		# ***case 1: power-law model***
#    	fit.res <- lm(z~dist+dd,data=td)		# ***case 2: exponential model***
        ssr <- sum(fit.res$residuals^2)
	return(ssr)
}


fi2 <- function(obs, sv, wind_dir){														
	td <- as.data.frame(sp_par_cal(obs, sv, wind_dir))
	colnames(td) <- c("x","y","dist","dd","cd", "dif","zr")
	td$z = log10(td$zr)
      	fit.res <- lm(z~dif+log(dist),data=td)		# ***case 1: power-law model ***
#       fit.res <- lm(z~dist+dd,data=td)		# ***case 2: exponential model***
     	ssr <- sum(fit.res$residuals^2)
	return(c(ssr, fit.res$coefficients, summary(fit.res)$r.squared))
}

cg <- function(sv, h){					
	output  = cbind(c(rep(sv[,1],4)+sin(0:3*90*pi/180)*h),
		    c(rep(sv[,2],4)+cos(0:3*90*pi/180)*h))
	return(output)
}

baf_pre <- function(sv, obs, nlps){
	pot_ang <- as.data.frame((1:36)*10)									
	pot_ssr <- apply(pot_ang, MARGIN = 1, FUN=fi, obs = obs, sv = sv)
	lowest_rank = which(rank(pot_ssr)==1)			# ***Power-law***
	next_lowest_rank = which(rank(pot_ssr)==2)		# ***Power-law***
	#lowest_rank = min(which(rank(pot_ssr)<2))				# ***exponential***
	#next_lowest_rank = min(which(rank(pot_ssr)<4 & rank(pot_ssr)>2))	# ***exponential***
	if(abs(lowest_rank - next_lowest_rank)==1 || abs(lowest_rank - next_lowest_rank)==35){
		return(c(1, pot_ang[lowest_rank,1]))
	}else{
		return(c(2, pot_ang[lowest_rank,1], pot_ang[next_lowest_rank, 1]))
	}
}

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

# Reference:
#	Gonzalez-Mellado, A. O., and S. De la Cruz-Reyna. "A simple semi-empirical approach to model thickness of 
# ash-deposits for different #eruption scenarios." Natural Hazards and Earth System Sciences 10.11 (2010): 2241.
# 	Yang Q, Bursik M. A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall 
# deposits. Bulletin of Volcanology, 2016, 78(10): 68.



