YMAT <- generateSim(1000)

ptm <- proc.time()

sq1<-SquareHinge(YMAT$YMAT)
#h1<-Hinge(YMAT$YMAT)
l1<-Logistic(YMAT$YMAT)

x1<-proc.time() - ptm

ptm <- proc.time()

sq2<-SquareHingeC(YMAT$YMAT,returnAll =T)
h2<-HingeC(YMAT$YMAT,returnAll =T)
l2<-LogisticC(YMAT$YMAT,returnAll =F)

x2<-proc.time() - ptm

as.numeric(x1)[1]/as.numeric(x2)[1]


library(animation)


## set some options first 


ani.options(interval = 0.1, nmax = floor(sq2$NN/block_size))
## use a loop to create images one by one
for (i in 1:ani.options('nmax')) {
  
  j=i*block_size;
  
  plot(YMAT$XX[1:j,], pch=(YMAT$YY+2)[1:j],col=heat.colors(sq2$NN)[1:j], ylim=c(min(YMAT$XX[,2]),max(YMAT$XX[,2])), xlim=c(min(YMAT$XX[,1]),max(YMAT$XX[,1])), xlab="X1",ylab="X2" )
  
  curve(-sq2$THETA_list[j,1]/sq2$THETA_list[j,3]-sq2$THETA_list[j,2]/sq2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
  #curve(-h2$THETA_list[j,1]/h2$THETA_list[j,3]-h2$THETA_list[j,2]/h2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
  curve(-l2$THETA_list[j,1]/l2$THETA_list[j,3]-l2$THETA_list[j,2]/l2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
  
  ani.pause()   ## pause for a while ('interval')
}

block_size = 100;
saveGIF({
  ani.options( nmax = floor(sq2$NN/block_size))
  for (i in 1:ani.options('nmax')) {
    j=i*block_size;
    
    plot(YMAT$XX[1:j,], pch=(YMAT$YY+2)[1:j],col=heat.colors(sq2$NN)[1:j], ylim=c(min(YMAT$XX[,2]),max(YMAT$XX[,2])), xlim=c(min(YMAT$XX[,1]),max(YMAT$XX[,1])), xlab="X1",ylab="X2" )
    
    curve(-sq2$THETA_list[j,1]/sq2$THETA_list[j,3]-sq2$THETA_list[j,2]/sq2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
    #curve(-h2$THETA_list[j,1]/h2$THETA_list[j,3]-h2$THETA_list[j,2]/h2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
    curve(-l2$THETA_list[j,1]/l2$THETA_list[j,3]-l2$THETA_list[j,2]/l2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
    }

}, interval = 0.1, movie.name = 'cool_demo.gif', ani.width = 600, ani.height = 600)



 plot(YMAT$XX,pch=(YMAT$YY+2),col=heat.colors(sq1$NN)[1:sq1$NN])
 curve(-sq1$THETA[1]/sq1$THETA[3]-sq1$THETA[2]/sq1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
 curve(-h1$THETA[1]/h1$THETA[3]-h1$THETA[2]/h1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
 curve(-l1$THETA[1]/l1$THETA[3]-l1$THETA[2]/l1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
 