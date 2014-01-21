ttestit=function(x,ncolumns,nreps,vbar,w)
  #assumes asin-sqrt transform before passing it here!
  {
    x1=x[4:28][1:nreps]
    x2=x[29:53][1:nreps]
    m1=mean(x1)
    m2=mean(x2)
    v1=var(x1)
    v2=var(x2)
    #corrected 7/31/2013
    t = (m1-m2)/(sqrt( ((1-w)/nreps)*(v1+v2) + 2*w*vbar/nreps ))
    return (t)
  }

x=read.table(gzfile(fn))
x[,4:ncol(x)]=asin(sqrt(x[,4:ncol(x)]))
NC=25
NREP=c(2,5,10,15,25)


results = matrix(data=NA,nrow=nrow(x),ncol=3+length(NREP),
  dimnames=list(c(rep(NA,times=nrow(x))),c("pos","ifreq","issel","r2","r5","r10","r15","r25")))

#predDF=read.table("/home/kevin/work/sim_expevol/parameter_space/predict_DF_using_anova/predictedDF.txt",header=T)

for( i in 1:length(NREP) )
  {
    D = c(rep(0,3),rep(1,NREP[i]),rep(0,NC-NREP[i]),rep(2,NREP[i]),rep(0,NC-NREP[i]))
    vbar = mean( c(apply( x[,D==1],1,var),apply( x[,D==2],1,var)) )
    x.tt = apply(x,1,ttestit,NC,NREP[i],vbar,0.1)
                                        #store the -log of the p-value
                                        #results[,i+3]=-pt(x.tt,df=(NREP[i]-1),log.p=TRUE,lower.tail=FALSE)
    #pDF = predDF$DFhat[which(predDF$r == NREP[i] & predDF$Nc == Nc)]
    results[,i+3]=x.tt#-(log(2)+pt(abs(x.tt),df=pDF,log.p=TRUE,lower.tail=FALSE))
  }
results[,1]=x[,1]
results[,2]=x[,2]
results[,3]=x[,3]
if(file.exists(ofn))
  {
    file.remove(ofn)
  }
write.table(results,file=ofn,row.names=FALSE,col.names=TRUE,quote=FALSE)
