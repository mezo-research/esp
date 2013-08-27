#####################################################################################################################
#
# 	espf.r
#	R code defing estimation and plot functions for use with implementation script (esp.r).
#
# 	Code by Athol Whitten (athol.whitten@mezo.com.au)
# 	Melbourne, Australia, 2012
#
#	Created: 16th May 2012
#	Updated: 10th November 2012
#
#######################################################################################################################

#Define Selectivity Function (Gamma): Relative selectivities with functional form of the standard Gamma function (with modal value rescaled to one):
esp.gamma <- function(length.class,gsize,theta){
	B=-0.5*(theta[1]*gsize-((theta[1]^2)*(gsize^2)+4*theta[2])^0.5)
	A=(theta[1]*gsize)/B
	sel=((length.class/(A*B))^A)*exp(A-(length.class/B))
	return(sel)
	}

#Define Selectivity Function (LNorm): Relative selectivities with functional form of the standard Log-Normal function (with modal value rescaled to one):	
esp.lnorm <- function(length.class,theta){
	sel=exp((-(log(length.class)-theta[3])^2)/(2*theta[4]^2) - (theta[4]^2)/2 + theta[3])/length.class
	return(sel)
	}


#Define a log-likelihood function to be maximised (implemented as a negative log-likelihood and minimised):
esp.nll <- function(theta,cdata,gsizes,sel.a,sel.b){
	catch.data=cdata[which(apply(cdata[,-1],1,sum)>0),-1]
	length.class=cdata[which(apply(cdata[,-1],1,sum)>0),1]
	sels.a=outer(length.class,gsizes,sel.a,theta)
	sels.b=sel.b(length.class,theta) # This part will probably have to be extended to use of the outer product too, which should would over single vectors?
	smatrix=cbind(sels.a,sels.b)
	mu=apply(catch.data,1,sum)/apply(smatrix,1,sum)
	esp.nll=-sum(catch.data*(log(mu*smatrix))-mu*smatrix)
	return(esp.nll)
	}	

#Create function to get estimates from particular fit:
esp.est <- function(fit,gsizes){
	theta=fit$par
	B1=-0.5*(theta[1]*gsizes[1]-((theta[1]^2)*(gsizes[1]^2)+4*theta[2])^0.5)
	A1=(theta[1]*gsizes[1])/B1
	B2=-0.5*(theta[1]*gsizes[2]-((theta[1]^2)*(gsizes[2]^2)+4*theta[2])^0.5)
	A2=(theta[1]*gsizes[2])/B2
	Net1.Mode=A1*B1
	Net2.Mode=A2*B2
	Net.Variance=(A1+1)*B1^2
	Net.StDev=sqrt(Net.Variance)
	b.Mode=exp(theta[3]-(theta[4]^2))
	b.Variance=exp(2*theta[3]+(theta[4]^2))*(exp(theta[4]^2)-1)
	parmest=rbind("Theta 1"=fit$par[1],"Theta 2"=fit$par[2],"Theta 3"=fit$par[3],"Theta 4"=fit$par[4])
	stnderr=sqrt(diag(solve(fit$hessian))) #The standard error of a parameter estimates is the square root of the diagonal of the inverse of the hessian matrix.
	output.a=cbind(parmest,stnderr)
	colnames(output.a)=c("Estimate","S.E.")
	output.b=rbind("Net Mode (Size 1)"=Net1.Mode,"Net Mode (Size 2)"=Net2.Mode,"Net Variance"=Net.Variance,"Net StDev"=Net.StDev,"B Mode"=b.Mode,"B Variance"=b.Variance)
	colnames(output.b)="Estimate"
	print(output.a)
	print(output.b)
	}

#Create function to get plots of data and of selectivity curves;
esp.plot <- function(fit,gsizes,cdata,sel.a,sel.b,plot.lens,label=TRUE,save=FALSE,save.to="directory",name="new",BS="TRUE"){
		
		theta=fit$par
		sels.a=outer(plot.lens,gsizes,sel.a,theta)
		sels.b=sel.b(plot.lens,theta)
		smatrix=cbind(plot.lens,sels.a,sels.b)
		
		windows(record=TRUE,width=900,height=500)
		
		barplot(as.matrix(cdata[,-1]),beside=TRUE,ylab="Frequency",xlab="Gear Type",ylim=c(0,1.1*max(cdata[,-1])))
		
		if(BS==TRUE){
			matplot(plot.lens,sels.a[,-1],type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
			matlines(plot.lens,sels.a[,1],type="l",col=1,lty=3)
			}
		
		if(BS==FALSE){
			matplot(plot.lens,sels.a,type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
			}
			
		matlines(plot.lens,sels.b,type="l",col=1,lwd=2)
			
			if(label==TRUE){
				for(i in 1:ncol(cdata[-1])){
				text(smatrix[which(smatrix[,i+1]==max(smatrix[,i+1])),1],1.03,colnames(cdata[-1])[i],cex=0.8)
				}
			}
		
		if(save==TRUE){
				png(file=paste(save.to,"/",name,".data.png",sep=""),width=900,height=500)
				barplot(as.matrix(cdata[,-1]),beside=TRUE,ylab="Frequency",xlab="Gear Type",ylim=c(0,1.1*max(cdata[,-1])))
				dev.off()
				
				png(file=paste(save.to,"/",name,".selectivity.png",sep=""),width=900,height=500)
				
				
					if(BS==TRUE){
						matplot(plot.lens,sels.a[,-1],type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
						matlines(plot.lens,sels.a[,1],type="l",col=1,lty=3)
						}
		
					if(BS==FALSE){
					matplot(plot.lens,sels.a,type="l",lty=2,col=1,xlab="Length (cm)", ylab="Relative Selectivity",ylim=c(0,1.05))
					}
				
				
				matlines(plot.lens,sels.b,type="l",col=1,lwd=2)
					
					if(label==TRUE){
						for(i in 1:ncol(cdata[-1])){
						text(smatrix[which(smatrix[,i+1]==max(smatrix[,i+1])),1],1.03,colnames(cdata[-1])[i],cex=0.8)
						}
					}
				dev.off()
				}
	}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This code/program has been developed by Mezo Research, Melbourne, Australia, 2012.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (http://www.fsf.org/).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++