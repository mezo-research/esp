#####################################################################################################################
#
# 	esp_perch.r
#	R code implementing selectivity estimation functions defined in espf.r
#	Plots and analyses of outputs for Perch Data.
#	
# 	Code by Athol Whitten (athol.whitten@mezo.com.au)
# 	Melbourne, Australia, 2011
#
#	Created: 16th May 2012
#	Updated: 19th November 2012
#
#######################################################################################################################

#Remove all objects and lists currently in workspace as precaution:
rm(list=ls())

#Specify directory of working folder;
main.folder <- "C:/ESP"

#Specify folder for printing plots of results (set here to be a subfolder of the working folder, and named 'Output'):
output.folder <- paste(main.folder,"/Output",sep="")
dir.create(output.folder)

#Source the requisite functions from the espf.r file:
source(paste(main.folder,"/Program/espf.r",sep=""))

#Get the data and set the gsizes (gear sizes) vector (perch used here as an example):
perch.data <- (read.table(paste(main.folder,"/Data/Perch.dat",sep=""),header=TRUE))
perch.gsizes <- c(5,14,17,21.5,25,30,33,38,45)

#View and check the data and gear size vector:
print(perch.data)
print(perch.gsizes)

#Set initial parameter values for Theta 1 through Theta 4;
gn.parms <- c(1,10) #Theta 1 and 2.
gc.parms <- c(2,0.5) #Theta 3 and 4.

#Implement optimisation with all data, using Nelder-Mead method as part of R optim function:
perch.fit <- optim(c(gn.parms,gc.parms),esp.nll,cdata=perch.data,gsizes=perch.gsizes,sel.a=esp.gamma,sel.b=esp.lnorm,hessian=TRUE,method="Nelder-Mead")

#Get parameter estimates and the mode and variance estimates for the first two nets and the Cormorants;
perch.est <- esp.est(fit=perch.fit,gsizes=perch.gsizes)

#Get plots of base data (barplot) and of the estimated selectivity curves over a specified range of plot lengths;
esp.plot(fit=perch.fit,gsizes=perch.gsizes,cdata=perch.data,sel.a=esp.gamma,sel.b=esp.lnorm,BS=TRUE,plot.lens=seq(0.0001,50,0.01),label=TRUE,save=TRUE,save.to=output.folder,name="perch")

#######################################################################################################################
# End of esp.r.
#######################################################################################################################