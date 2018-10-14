#-------------------------------------------------#
# FIGURE
# Association Between Enactment of a “Stand Your 
# Ground” Self-defense Law and Unlawful Homicides 
# in Florida
# JAMA Internal Medicine
# DOI:10.1001/jamainternmed.2017.3433
# 
# David K. Humphreys, Antonio Gasparrini, 
# Douglas J. Wiebe (2017b)
# Email: david.humphreys@spi.ox.ac.uk
#-------------------------------------------------#

tiff(file="JAMA_IM_Figure",width=4000,height=3000,res=600)
par (mar=c(4,4,1.5,1.5))
jhom.obs <- with(jhom, unlaw.h/stdpop*10^5)
jhom.datanew <- data.frame(stdpop=mean(jhom$stdpop),Effective=rep(c(0,1),c(819,1221)),
                             time= 1:2040/10,month=rep(1:120/10,17))
jlawhom.m3.pred1 <- predict(jlawhom.m3,type="response",jhom.datanew)/mean(jhom$stdpop)*10^5
jlawhom.m3.pred2 <- predict(jlawhom.m3,type="response",transform(jhom.datanew,month=4.5))/mean(jhom$stdpop)*10^5

plot(1:204,jhom.obs,type="n",ylim=c(0.0, 1.0),xlab="Year",
     ylab="Rate per 100,000",frame.plot=F,xaxt="n",las=2) #specifies the plot area
rect(82, 0.0, 204, 1.0, col=grey(0.9),border=F) # format of the rectangle (xleft, ybottom, xright, ytop)
points(1:204,jhom.obs,cex=0.7, pch=16, col="dodgerblue4") # Points 1:59 of obs - but there are only 59 points
axis(1,at=0:17*12,labels=F) #The rule of x axis- no labes/ 5 sections
axis(1,at=0:16*12+6,tick=F,labels=1999:2015)  # What goes in them
lines(1:2040/10,jlawhom.m3.pred1,col="darkgoldenrod3", lwd=2)
lines(1:2040/10,jlawhom.m3.pred2,col="darkgoldenrod3",lty=2)
title("Unlawful Homicide in Florida 1999-2015")
legend ("topleft", legend=c("Linear effect", "Seasonally adjusted effect"), col=c("darkgoldenrod3", "darkgoldenrod3"), lty=2:1, cex=0.7, bty="n", inset=0.015)
dev.off()
