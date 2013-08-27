#####################################################################################################################
#
# 	esp.r
#	R code implementing selectivity estimation functions defined in espf.r
#	Plots and analyses of outputs provided in this script.
#	
# 	Code by Athol Whitten (athol.whitten@mezo.com.au)
# 	Melbourne, Australia, 2012
#
#	Created: 16th May 2012
#	Updated: 19th November 2012
#
#######################################################################################################################

#Remove all objects and lists currently in workspace as precaution:
rm(list=ls())

#Specify directory of working folder;
main.folder <- "C:/Dropbox/Mezo/Shared/Projects/Esp"

#Specify folder for printing plots of results (set here to be a subfolder of the working folder, and named 'Output'):
output.folder <- paste(main.folder,"/Output",sep="")
dir.create(output.folder)

#Source the requisite functions from the espf.r file:
source(paste(main.folder,"/Program/espf.r",sep=""))

#Get the data and set the gsizes (gear sizes) vector:
species.data <- (read.table(paste(main.folder,"/Data/Species.dat",sep=""),header=TRUE))
species.gsizes <- c(1,2,3,4,5,6,7)

#View and check the data and gear size vector:
print(species.data)
print(species.gsizes)

#Set initial parameter values for Theta 1 through Theta 4;
net.parms <- c(1,2) #Theta 1 and 2.
b.parms <- c(3,4) #Theta 3 and 4.

#Implement optimisation with all data, using Nelder-Mead method as part of R optim function:
species.fit <- optim(c(net.parms,b.parms),esp.nll,cdata=species.data,gsizes=species.gsizes,sel.a=esp.gamma,sel.b=esp.lnorm,hessian=TRUE,method="Nelder-Mead")

#Get parameter estimates and the mode and variance estimates for the first two nets and the Cormorants;
species.est <- esp.est(fit=species.fit,gsizes=species.gsizes)

#Get plots of base data (barplot) and of the estimated selectivity curves over a specified range of plot lengths;
esp.plot(fit=species.fit,gsizes=species.gsizes,cdata=species.data,sel.a=esp.gamma,sel.b=esp.lnorm,BS=TRUE,plot.lens=seq(0.0001,50,0.01),label=TRUE,save=TRUE,save.to=output.folder,name="species")

#######################################################################################################################
# End of esp.r.
#######################################################################################################################