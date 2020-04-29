check.integer <- function(x) {
  x == round(x)
}
SEIR.model<- function(current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S  <- state_values [1]        # susceptibles
  E  <- state_values [2]        # exposed
  I  <- state_values [3]        # infectious
  R  <- state_values [4]        # recovered
  input <- approxfun(parameters$beta, rule = 2)
  with(as.list(c(state_values, parameters)), {
    betat <-input(current_timepoint)
    dS = -(betat * S * I/N)-mu*(1-S)
    dE = (betat * S * I/N) - (sigma+mu) * E
    dI = (sigma * E) - (gamma+mu) * I
    dR = (gamma * I)-mu*R
    results = c(dS, dE, dI, dR)
    return(list(results))
  })
}

SEIR_TV<-function(data,n,par_start=NULL,theta0,status0,nknots,ahead,loss="poisson"){
#define splines BASIS
X.new<-ns(0:(n+ahead),nknots,Boundary.knots = range(0:(n+7)),intercept = TRUE)
X<-X.new[1:n,]
obs<-data$nuovi_positivi
obs2<-data$totale_positivi
#its important to specify a correct initial starting point
if(is.null(par_start)) par_start<-c(-3,seq(0,-4,len=ncol(X)))
est<-nlm(f=f1,p=par_start,gradtol = 1e-7, 
         parameters=theta0 ,initials=status0,
         obs=obs,X=X,type=loss,obs2=obs2)  
return(list(fit=est,X=X,X.new=X.new,obs=obs,obs2=obs2))
}

f1_vec<-function(par,parameters ,initials,obs,X, type="wls",obs2=NULL){
  gammap<-inv.logit(par[1])
  betat<-inv.logit(as.numeric(X%*%par[-1]))
  theta.list<-list(mu =parameters[1], sigma = parameters[2], gamma =gammap,N=parameters[4],beta=betat)
    n<-length(obs) 
  times <- seq(0, n, by=0.1) 
  out <- lsode(y=initials, times =times,func= SEIR.model, parms =theta.list)
    expected<-out[check.integer(out[,1]),"E"]*parameters[2]
    expected2<-out[check.integer(out[,1]),"I"]
    if (type == "ls")
      v<-((obs-expected[-n])^2+(obs2-expected2[-1])^2)
    if (type == "wls")
      v<-(((obs-expected[-n])^2)/expected[-n]+((obs2-expected2[-1])^2)/expected2[-1])
    if (type == "poisson")
      v<- (-dpois(obs,expected[-n],log = TRUE) -dpois(obs2,expected2[-1],log = TRUE))
    return(v)
}


f1<-function(par,parameters,initials,obs,X,type="wls",obs2)
{
  v<-f1_vec(par,parameters=parameters ,initials=initials,
            obs=obs,X=X,type=type,obs2=obs2)
  sum(v)
}




HAC.weights<-function(n,lag=0,kernel.type="NW")
{
  if (kernel.type == "NW")
  {
    l<-0:(n-1)
    w<-(1-l/(lag+1))*(l <= (lag+1))
  }
  if (kernel.type == "Andrews"){
    w<-numeric(n)
    w[1]<-1
    l<-(6*pi/5)*(1:(n-1))/(lag+1)
    w[-1]<-(3/l^2)*(sin(l)/l-cos(l))
  }
  return(w)      
}
sandwich<-function(fitted.model, parameters,initials,obs,X, type ="wls",lag=1,kernel.type="NW",obs2)
{ 
  if (!is.null(obs))
    n<-length(obs) 
  else
    n<-length(obs2) 
  
  
  x<-fitted.model$estimate
  p<-length(x)
  g<-jacobian(func=f1_vec,x=x, parameters=parameters ,
              initials=initials,obs=obs,X=X,type=type,obs2=obs2)
  g<-as.matrix(g)
  # Hn
  Hn<-numDeriv::hessian(func=f1,x=x, 
                        parameters=parameters ,initials=initials,
                        obs=obs,X=X,type=type,obs2=obs2)
  Hn<-Hn/n
  # Jn
  Jn<-crossprod(g)/n
  w<-HAC.weights(n,lag = lag,kernel.type = kernel.type)
  Jn.hac <- 0.5 * crossprod(g) * w[1] 
  for(i in 2:length(w)) 
    Jn.hac <- Jn.hac + w[i]*t(matrix(g[1:(n-i+1),],ncol=p))%*%g[i:n,]
  Jn.hac<-(Jn.hac + t(Jn.hac))/n
  Hn.inv<-solve(Hn)
  clic <- fitted.model$minimum - sum(diag(Hn.inv%*%Jn))
  clic.hac <- fitted.model$minimum - sum(diag(Hn.inv%*%Jn))
  v<-Hn.inv%*%Jn%*%t(Hn.inv)/n
  v.hac<-Hn.inv%*%Jn.hac%*%t(Hn.inv)/n
  return(list(est=x,g=g,Hn=Hn,Jn=Jn,varcov=v,clic=clic,Jn.hac=Jn.hac,varcov.hac=v.hac,clic.hac=clic.hac))
}



