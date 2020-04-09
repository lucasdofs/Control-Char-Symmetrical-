############### Author: LUCAS SALES ############
#Institution: Universidade Federal do Rio Grande do Norte
#e-mail: lucasofsales@gmail.com  or ldo-sales@hotmail.com

#Symmetrical Class simulation #######

### This code was developed for distribution which belongs to symmetrical
#clas. More especifically, the code was illustrated for the t-Student and
#power-exponential distributions.
#However, to use other distributions of this class, change the rlnorm 
#to the desired distribution.






############  t distribution ########
#FUNCTIONS 
#m is the sample size
#phi the dispersion parameter
#mu the mean of the distribution
#nu is the degrees of freedom of the t-Student
#CL is an vector that contains the control limits CL=c(LCL,UCL)
# Bootstrap Control Limits
t_limits=function(n,mu,phi,nu){
  UCL_t=LCL_t=NULL
  for(i in 1:10000){
    print(i)
    x=matrix(rt(5000*n,df=nu)*sqrt(phi)+mu,ncol=n)
    x=rowMeans(x)
    UCL_t[i]=quantile(x,1-0.0027/2)
    LCL_t[i]=quantile(x,0.0027/2)
  }
  (UCL_t=mean(UCL_t))
  (LCL_t=mean(LCL_t))
  return(c(LCL_t,UCL_t))
}

#Providing the estimated ARL 
t_arl=function(CL,n,mu,phi,nu){
  LCL_t=round(CL[1],2);UCL_t=round(CL[2],2)
  xi=nu/(nu-2)
  normal_CL=mu+3*c(-1,1)*sqrt(phi*xi/n) #limite de controle normal
  nmaf_n=nmaf_t=nma1_t=nma2_t=nma3_t=nma4_t=nma5_t=nma6_t=NULL
  for(i in 1:10000){
    print(i)
    x=matrix(rt(5000*n,df=nu)*sqrt(phi)+mu,ncol=n)
    x=rowMeans(x)
    nmaf_t[i]=which(x>UCL_t| x<LCL_t)[1]
    nmaf_n[i]=which(x>normal_CL[2]| x<normal_CL[1])[1]
    #x1=x-1.5* sqrt(phi*xi) #Shifted sample
    #nma1_t[i]=which(x1>UCL_t| x1<LCL_t)[1]
    #x2=x-2*sqrt(phi*xi)    #Shifted sample
    #nma2_t[i]=which(x2>UCL_t| x2<LCL_t)[1]
    #x3=x-3*sqrt(phi*xi)    #Shifted sample
    #nma3_t[i]=which(x3>UCL_t| x3<LCL_t)[1]
    #x4=x+1.5* sqrt(phi*xi) #Shifted sample
    #nma4_t[i]=which(x4>UCL_t| x4<LCL_t)[1]
    #x5=x+2*sqrt(phi*xi)    #Shifted sample
    #nma5_t[i]=which(x5>UCL_t| x5<LCL_t)[1]
    #x6=x+3*sqrt(phi*xi)    #Shifted sample
    #nma6_t[i]=which(x6>UCL_t| x6<LCL_t)[1]
  }
  aux=data.frame("ARL0_Normal"=mean(nmaf_n,na.rm=T),"ARL0_T"=mean(nmaf_t,na.rm=T),
                 "delta= -1.5"=mean(nma1_t),"delta= -2"=mean(nma2_t),
                 "delta= -3"=mean(nma3_t),"delta= 1.5"=mean(nma4_t),"delta= 2"=mean(nma5_t),
                 "delta= 3"=mean(nma6_t))
  return(aux)
}


############## Power Exponential 
#FUNCTIONS 
#m is the sample size
#phi the dispersion parameter
#mu the mean of the distribution
#k is the parameters of shape ( -1<k<1)
#CL is an vector that contains the control limits CL=c(LCL,UCL)


#Generation of a power-exponential sample
rpower.exp=function(m,mu,phi,k){
  r=2/(1+k)
  v=runif(m,-1,1)
  w=rgamma(m,1+1/r,1)
  z= (2*w)^(1/r)*v
  y= z*sqrt(phi)+mu
  return(y)
}

# Bootstrap Control Limits
PE_limits=function(n,mu,phi,k){
  UCL_PE=LCL_PE=NULL
  for(i in 1:10000){
    print(i)
    x=matrix(rpower.exp(5000*n,mu,phi,k),ncol=n)
    x=rowMeans(x)
    UCL_PE[i]=quantile(x,1-0.0027/2)
    LCL_PE[i]=quantile(x,0.0027/2)
  }
  (UCL_PE=mean(UCL_PE))
  (LCL_PE=mean(LCL_PE))
  return(c(LCL_PE,UCL_PE))
}

#Estimated ARL
PE_arl=function(CL,n,mu,phi,k){
  LCL_PE=CL[1];UCL_PE=CL[2]
  xi=2^(1+k)* gamma(1.5*(1+k))/ gamma((1+k)/2)
  normal_cl=mu+3*c(-1,1)*sqrt(phi*xi/n) #limite de controle normal
  nmaf_n=nmaf_PE=nma1_PE=nma2_PE=nma3_PE=nma4_PE=nma5_PE=nma6_PE=NULL
  nma1_n=nma2_n=nma3_n=nma4_n=nma5_n=nma6_n=NULL
  for(i in 1:10000){
    print(i)
    x=matrix(rpower.exp(5000*n,mu,phi,k),ncol=n)
    x=rowMeans(x)
    # normal_CL=mean(x)+3*c(-1,1)*sd(x)/sqrt(n)#limite de controle normal
    nmaf_PE[i]=which(x>UCL_PE| x<LCL_PE)[1]
    nmaf_n[i]=which(x>normal_cl[2]| x<normal_cl[1])[1]
    x1=x-1.5* sqrt(phi*xi) #Amostral alterada
    nma1_PE[i]=which(x1>UCL_PE| x1<LCL_PE)[1]
    nma1_n[i]=which(x1>normal_cl[2]| x1<normal_cl[1])[1]
    x2=x-2*sqrt(phi*xi)    #Amostral alterada
    nma2_PE[i]=which(x2>UCL_PE| x2<LCL_PE)[1]
    nma2_n[i]=which(x2>normal_cl[2]| x2<normal_cl[1])[1]
    x3=x-3*sqrt(phi*xi)    #Amostral alterada
    nma3_PE[i]=which(x3>UCL_PE| x3<LCL_PE)[1]
    nma3_n[i]=which(x3>normal_cl[2]| x3<normal_cl[1])[1]
    x4=x+1.5* sqrt(phi*xi) #Amostral alterada
    nma4_PE[i]=which(x4>UCL_PE| x4<LCL_PE)[1]
    nma4_n[i]=which(x4>normal_cl[2]| x4<normal_cl[1])[1]
    x5=x+2*sqrt(phi*xi)    #Amostral alterada
    nma5_PE[i]=which(x5>UCL_PE| x5<LCL_PE)[1]
    nma5_n[i]=which(x5>normal_cl[2]| x5<normal_cl[1])[1]
    x6=x+3*sqrt(phi*xi)    #Amostral alterada
    nma6_PE[i]=which(x6>UCL_PE| x6<LCL_PE)[1]
    nma6_n[i]=which(x6>normal_cl[2]| x6<normal_cl[1])[1]
    
  }
  aux=data.frame("ARL0_Normal"=mean(nmaf_n,na.rm=T),"ARL0_PE"=mean(nmaf_PE,na.rm=T),
                 "delta= -1.5"=mean(nma1_PE),"delta= -2"=mean(nma2_PE),
                 "delta= -3"=mean(nma3_PE),"delta= 1.5"=mean(nma4_PE),"delta= 2"=mean(nma5_PE),
                 "delta= 3"=mean(nma6_PE),
                 "delta= -1.5"=mean(nma1_n),"delta= -2"=mean(nma2_n),
                 "delta= -3"=mean(nma3_n),"delta= 1.5"=mean(nma4_n),"delta= 2"=mean(nma5_n),
                 "delta= 3"=mean(nma6_n))
  return(aux)
}
