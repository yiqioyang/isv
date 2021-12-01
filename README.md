This tool presents a new method to identify the vent location of tephra fall deposits based on thickness or maximum clast size measurements. It is temporarily named "svl", which is short for Source Vent Locator. The method estimates the vent location by coupling semi-empirical models of tephra thickness distribution with a gradient descent method. Two distinct semi-empirical models, proposed by Gonzalez-Mellado and De la Cruz-Reyna (2010) and Yang and Bursik (2016), can be used for the method. Based on their characteristics, they are named the power-law and exponential models, respectively. The two semi-empirical models assume different thinning rate with distance to the source vent.

This offline tool includes three R scripts, and thickness dataset of North Mono Bed 2 digitized from Sieh and Bursik (1986).
The files are:
	pub_svl_2.0_exponential.r:
		Functions of the method "svl" that are coupled with the semi-empirical model proposed by Yang and Bursik (2016). 

	pub_svl_2.0_power-law.r
		Functions of the method "svl" that are coupled with the semi-empirical model proposed by Gonzalez-Mellado and De la Cruz-Reyna (2010).
		This script is only slightly different from "pub_svl_2.0_exponential.r". The differences are marked with "***" on the right side of the corresponding lines.
	
	demo_svl_2.0.r
		Demonstration of how this method works, using data for North Mono Bed 2. 		

Our experiments show that the method is able to give an accurate estimate of the source vent location of tephra fall deposits with > 10 evenly-distributed sample sites. The accuracy is affected by characteristics of the tephra deposit, quality of input dataset, and features of the semi-empirical model being adopted. With >20 input points, the semi-empirical model of Gonzalez-Mellado and De la Cruz-Reyna (2010) provides a more accurate estimate. With 6-20 input points, we recommend the use of the semi-empirical model of Yang and Bursik (2016), owing to its stability. 

Our experiments also indicate that in the case of sparse observations (6-~15), it is more robust to constrain the vent location by estimating the dispersal axis (a line that is defined by the estimated source vent and prevailing wind direction) using the semi-empirical model of Yang and Bursik (2016). This is because estimating the dispersal axis only relies on the assumption that the wind direction is constant, and is not affected by how the deposit thins with distance from the vent. 

The manuscript introducing this method has been submitted to the Bulletin of Volcanology. 
We appreciate your comments, suggestions, and feedback. We hope that our method could be helpful to your work. 

*** For users who are interested in trying the method with the example dataset, or their own dataset, we recommend them to read the notations in "demo_svl_2.0.r". The notations in it include a self-consistent and brief description of how the method works, and how to implement the method. For users who are interested in how the method works in greater detail, please refer to the notes and scripts in "pub_svl_2.0_exponential.r".



References:
	Gonzalez-Mellado A O, De la Cruz-Reyna S. A simple semi-empirical approach to model thickness of ash-deposits for different eruption scenarios. Natural Hazards and Earth System Sciences, 2010, 10(11): 2241.
	Sieh, Kerry, and Marcus Bursik. "Most recent eruption of the Mono Craters, eastern central California." Journal of Geophysical Research: Solid Earth 91.B12 (1986): 12539-12571.
	Yang Q, Bursik M. A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall deposits. Bulletin of Volcanology, 2016, 78(10): 68.
