library(SSOSVM)
library(caret)
library(animation)
###########################################

YMAT <- generateSim(10^4)

ptm <- proc.time()

sq1<-SquareHinge(YMAT$YMAT,returnAll =T)
h1<-Hinge(YMAT$YMAT,returnAll =T)
l1<-Logistic(YMAT$YMAT,returnAll =T)

x1<-proc.time() - ptm

ptm <- proc.time()
sq2<-SquareHingeC(YMAT$YMAT,returnAll =T)
h2<-HingeC(YMAT$YMAT,returnAll =T)
l2<-LogisticC(YMAT$YMAT,returnAll =T)
x2<-proc.time() - ptm
x2

as.numeric(x1)[1]/as.numeric(x2)[1]
############################################

saveGIF({
  block_size = floor(sq2$NN/100);
  ani.options( nmax = floor(sq2$NN/block_size))
  ## use a loop to create images one by one
  for (i in 1:ani.options('nmax')) {
    
    j=i*block_size;
    
    plot(YMAT$XX[1:j,], pch=(YMAT$YY+2)[1:j],col=heat.colors(sq2$NN)[1:j], ylim=c(min(YMAT$XX[,2]),max(YMAT$XX[,2])), xlim=c(min(YMAT$XX[,1]),max(YMAT$XX[,1])), xlab="X1",ylab="X2" )
    
    curve(-sq2$THETA_list[j,1]/sq2$THETA_list[j,3]-sq2$THETA_list[j,2]/sq2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
    curve(-h2$THETA_list[j,1]/h2$THETA_list[j,3]-h2$THETA_list[j,2]/h2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
    curve(-l2$THETA_list[j,1]/l2$THETA_list[j,3]-l2$THETA_list[j,2]/l2$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
    
    ani.pause()   ## pause for a while ('interval')
  }
}, interval = 0.1, movie.name = 'cool.gif', ani.width = 300, ani.height = 300)


saveGIF({
  block_size = floor(sq2$NN/100);
  ani.options( nmax = floor(sq2$NN/block_size))
  ## use a loop to create images one by one
  for (i in 1:ani.options('nmax')) {
    
    j=i*block_size;
    
    plot(YMAT$XX[1:j,], pch=(YMAT$YY+2)[1:j],col=heat.colors(sq2$NN)[1:j], ylim=c(min(YMAT$XX[,2]),max(YMAT$XX[,2])), xlim=c(min(YMAT$XX[,1]),max(YMAT$XX[,1])), xlab="X1",ylab="X2" )
    
    curve(-sq1$THETA_list[j,1]/sq1$THETA_list[j,3]-sq1$THETA_list[j,2]/sq1$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
    curve(-h1$THETA_list[j,1]/h1$THETA_list[j,3]-h1$THETA_list[j,2]/h1$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
    curve(-l1$THETA_list[j,1]/l1$THETA_list[j,3]-l1$THETA_list[j,2]/l1$THETA_list[j,3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
    
    ani.pause()   ## pause for a while ('interval')
  }
}, interval = 0.1, movie.name = 'bad.gif', ani.width = 300, ani.height = 300)


############################################
# # download data from http://yann.lecun.com/exdb/mnist/
# download.file("http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz",
#               "train-images-idx3-ubyte.gz")
# download.file("http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz",
#               "train-labels-idx1-ubyte.gz")
# download.file("http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",
#               "t10k-images-idx3-ubyte.gz")
# download.file("http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz",
#               "t10k-labels-idx1-ubyte.gz")
# 
# # gunzip the files
# R.utils::gunzip("train-images-idx3-ubyte.gz")
# R.utils::gunzip("train-labels-idx1-ubyte.gz")
# R.utils::gunzip("t10k-images-idx3-ubyte.gz")
# R.utils::gunzip("t10k-labels-idx1-ubyte.gz")

# load images
train = load_image_file("train-images-idx3-ubyte")
test  = load_image_file("t10k-images-idx3-ubyte")

testlab  = (load_label_file("t10k-labels-idx1-ubyte"))
trainlab = (load_label_file("train-labels-idx1-ubyte"))

k=1
YY<-2.0*((trainlab==k)-.5)
YY2<-2.0*((testlab==k)-.5)



pca1<-princomp(train)
B=2
test1<-predict(pca1, newdata = test)[,1:B]
#set<-sample(c(which(trainlab==1), sample(which(trainlab!=1), length(which(trainlab==1)))
trainpc<-pca1$scores[set,1:B]
YMAT <- cbind(YY[set],as.matrix(trainpc))



M=20
set<-sample(60000,M)
B=784
YMAT <- cbind(YY[set],as.matrix(train[set,]))

ptm <- proc.time()
h2<-HingeC(YMAT, DIM = B, returnAll =F)
x2<-proc.time() - ptm
x2



xl<- sign(rowSums(t(apply(test1, 1 , `*` , l2$THETA[-(1)] )))+l2$THETA[(1)]) 
xsq<- sign(rowSums(t(apply(test1, 1 , `*` , sq2$THETA[-(1)] )))+sq2$THETA[(1)]) 
xh<- sign(rowSums(t(apply(test1, 1 , `*` , h2$THETA[-(1)] ))) +h2$THETA[(1)]) 

if(B==2){
  plot(trainpc, pch=(YY[set]+2),col=as.factor(YY[set]+2), ylim=c(min(trainpc[,2]),max(trainpc[,2])), xlim=c(min(trainpc[,1]),max(trainpc[,1])), xlab="X1",ylab="X2")
  curve(-sq2$THETA[1]/sum(sq2$THETA[-c(1,2)])-sq2$THETA[2]/sum(sq2$THETA[-c(1,2)])*x, col='orange',add=T)
  curve(-h2$THETA[1]/sum(h2$THETA[-c(1,2)])-h2$THETA[2]/sum(h2$THETA[-c(1,2)])*x,col='pink',add=T)
  curve(-l2$THETA[1]/sum(l2$THETA[-c(1,2)])-l2$THETA[2]/sum(l2$THETA[-c(1,2)])*x,col='green',add=T)
}

pred<-xl
truth<-YY2

xtab1 <- table(pred, truth)
confusionMatrix(xtab1)

pred<-xsq
xtab2<- table(pred, truth)
confusionMatrix(xtab2)


pred<-xh
xtab3<- table(pred, truth)
confusionMatrix(xtab3)





















## set some options first 




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
 curve(-sq2$THETA[1]/sq2$THETA[3]-sq2$THETA[2]/sq2$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
 
 curve(-h1$THETA[1]/h1$THETA[3]-h1$THETA[2]/h1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
 curve(-l1$THETA[1]/l1$THETA[3]-l1$THETA[2]/l1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
 