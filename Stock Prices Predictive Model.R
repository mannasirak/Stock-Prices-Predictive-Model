##Libraries##
library(fGarch)
library(sn)
library(tseries)
library(moments)
library(MASS)
library(quantmod)
library(copula)
library(LambertW)
library(ks)
library(evir)
library(forecast)
library(fracdiff)
library(cvar)


##Download data##
getSymbols(Symbols = 'MCD', src = 'yahoo', from = '2009-11-01', to = '2019-11-01')
getSymbols(Symbols = 'SPXL', src = 'yahoo', from = '2009-11-01', to = '2019-11-01')

##Creating and plotting time series of stock prices##
MCDprice=MCD[,6]
plot(MCDprice, main="MCD Stock Price", type='l')
summary(MCDprice)
SPXLprice=SPXL[,6]
plot(SPXLprice, main="SPXL Stock Price", type='l')
summary(SPXLprice)

##Creating log return variables##
MCDlogpr <- log(MCDprice)
MCDlogret <- diff(MCDlogpr)[-1,]
MCD <- MCDlogret
SPXLlogpr <- log(SPXLprice)
SPXLlogret <- diff(SPXLlogpr)[-1,]
SPXL <- SPXLlogret

##Exploring data## 
summary(MCD)
plot(MCD, main="MCD Stock Log Return", type="l")
lines(MCD)
summary(SPXL)
plot(SPXL, main="SPXL Stock Log Return", type="l")
lines(SPXL)

##Plotting log return ACFs##
acf(MCD, main="MCD Log Return ACF")
Box.test(x=MCD, lag=10, type="Ljung-Box")
acf(SPXL, main="SPXL Log Return ACF")
Box.test(x=SPXL, lag=10, type="Ljung-Box")

##Log return QQ plots##
qqnorm(MCD, main="MCD Log Return Normal QQ Plot")
qqline(MCD)
qqnorm(SPXL, main="SPXL Log Return Normal QQ Plot")
qqline(SPXL)

##Log return densities##
plot(density(MCD), main="MCD Density Stock Log Return", type="l")
plot(density(SPXL), main="SPXL Density Stock Log Return", type="l")

##Skewness and Kurtosis##
skewness(MCD)
skewness(SPXL)
kurtosis(MCD)
kurtosis(SPXL)

##Fitting skewed t-dist##
fit1=sstdFit(MCD,hessian=T)
fit1
fit2=sstdFit(SPXL,hessian=T)
fit2

##Fitting AR(1)-GARCH(1,1) model, normal error distribution##
adf.test(MCD)
adf.test(SPXL)

fitMCD=garchFit(formula=~arma(1,0)+garch(1,1),data=MCD,cond.dist='norm')
summary(fitMCD)
fitSPXL=garchFit(formula=~arma(1,0)+garch(1,1),data=SPXL,cond.dist='norm')
summary(fitSPXL)

##AIC and BIC##
aicMCD=2*(fitMCD@fit$llh)+2*length(fitMCD@fit$params$index)
aicSPXL=2*(fitSPXL@fit$llh)+2*length(fitSPXL@fit$params$index)
bicMCD=2*(fitMCD@fit$llh)+log(length(MCD))*length(fitMCD@fit$params$index)
bicSPXL=2*(fitSPXL@fit$llh)+log(length(SPXL))*length(fitSPXL@fit$params$index)
aicMCD
bicMCD
aicSPXL
bicSPXL

##Computing residuals##
MCDr=residuals(fitMCD)
MCDrst=residuals(fitMCD, standardize=TRUE)
acf(MCDr, main="MCD AR/GARCH Residuals ACF")
acf(MCDr^2, main="MCD AR/GARCH Squared Residuals ACF")
acf(MCDrst, main="MCD Standardized AR/GARCH Residuals ACF")
acf(MCDrst^2, main="MCD Standardized AR/GARCH Squared Residuals ACF")

SPXLr=residuals(fitSPXL)
SPXLrst=residuals(fitSPXL, standardize=TRUE)
acf(SPXLr, main="SPXL AR/GARCH Residuals ACF")
acf(SPXLr^2, main="SPXL AR/GARCH Squared Residuals ACF")
acf(SPXLrst, main="SPXL Standardized AR/GARCH Residuals ACF")
acf(SPXLrst^2, main="SPXL Standardized AR/GARCH Squared Residuals ACF")

##Testing for normality of standardized residuals##
qqnorm(MCDrst, main="MCD Normal QQ-Plot Stand. Residuals")
qqline(MCDrst)
jarque.bera.test(MCDrst)

qqnorm(SPXLrst, main="SPXL Normal QQ-Plot Stand. Residuals")
qqline(SPXLrst)
jarque.bera.test(SPXLrst)

##Fitting AR(1)-GARCH(1,1) model, t error distribution##
fitMCDtd=garchFit(formula=~arma(1,0)+garch(1,1),data=MCD,cond.dist='std')
summary(fitMCDtd)
fitSPXLtd=garchFit(formula=~arma(1,0)+garch(1,1),data=SPXL,cond.dist='std')
summary(fitSPXLtd)

##AIC and BIC, t error distribution##
aicMCD1=2*(fitMCDtd@fit$llh)+2*length(fitMCDtd@fit$params$index)
aicSPXL1=2*(fitSPXLtd@fit$llh)+2*length(fitSPXLtd@fit$params$index)
bicMCD1=2*(fitMCDtd@fit$llh)+log(length(MCD))*length(fitMCDtd@fit$params$index)
bicSPXL1=2*(fitSPXLtd@fit$llh)+log(length(SPXL))*length(fitSPXLtd@fit$params$index)
aicMCD1
bicMCD1
aicSPXL1
bicSPXL1

##Computing residuals and testing for t-dist, t error distribution##
MCDrtd=residuals(fitMCDtd)
MCDrsttd=residuals(fitMCDtd, standardize=TRUE)
acf(MCDrtd, main="MCD AR/GARCH Residuals ACF, t dist")
acf(MCDrtd^2, main="MCD AR/GARCH Squared Residuals ACF, t dist")
acf(MCDrsttd, main="MCD Standardized AR/GARCH Residuals ACF, t dist")
acf(MCDrsttd^2, main="MCD Standardized AR/GARCH Squared Residuals ACF, t dist")
ks.test.t(MCDrsttd)

SPXLrtd=residuals(fitSPXLtd)
SPXLrsttd=residuals(fitSPXLtd, standardize=TRUE)
acf(SPXLrtd, main="SPXL AR/GARCH Residuals ACF, t dist")
acf(SPXLrtd^2, main="SPXL AR/GARCH Squared Residuals ACF, t dist")
acf(SPXLrsttd, main="SPXL Standardized AR/GARCH Residuals ACF, t dist")
acf(SPXLrsttd^2, main="SPXL Standardized AR/GARCH Squared Residuals ACF, t dist")
ks.test.t(SPXLrsttd)

##Fitting parametric distribution to data##
plot(density(MCDrsttd), main="MCD Density Distribution With Standardized Residuals")
fitResidMCDtd=fitdistr(MCDrsttd, "t")
fitResidMCDsstd=sstdFit(MCDrsttd)
fitResidMCDtd$estimate
fitResidMCDsstd$estimate

plot(density(SPXLrsttd), main="SPXL Density Distribution With Standardized Residuals")
fitResidSPXLtd=fitdistr(SPXLrsttd, "t")
fitResidSPXLsstd=sstdFit(SPXLrsttd)
fitResidSPXLtd$estimate
fitResidSPXLsstd$estimate

##AIC and BIC, t error distribution (2)##
aicMCDtd=2*(fitResidMCDtd$loglik)+2*length(fitResidMCDtd$estimate)
bicMCDtd=2*(fitResidMCDtd$loglik)+log(length(MCD))*length(fitResidMCDtd$estimate)
aicMCDsstd=2*(fitResidMCDsstd$loglik)+2*length(fitResidMCDsstd$estimate)
bicMCDsstd=

aicSPXLtd=2*(fitResidSPXLtd$loglik)+2*length(fitResidSPXLtd$estimate)
bicSPXLtd=2*(fitResidSPXLtd$loglik)+log(length(SPXL))*length(fitResidSPXLtd$estimate)
aicSPXLsstd=2*(fitResidSPXLsstd$loglik)+2*length(fitResidSPXLsstd$estimate)
bicSPXLsstd=

aicMCDtd
aicMCDsstd
bicMCDtd
aicSPXLtd
aicSPXLsstd
bicSPXLtd

##Plots of t distributed error standardized residuals##
normResidMCD=(MCDrsttd-fitResidMCDtd$estimate[1])/fitResidMCDtd$estimate[2]
nMCD=length(normResidMCD)
x=qt((1:nMCD)/(nMCD+1),df=fitResidMCDtd$estimate[3])
qqplot(x,sort(normResidMCD),xlab="t",ylab="residuals", main = "MCD t QQ Plot Standardized Residuals")
qqline(normResidMCD)

normResidSPXL=(SPXLrsttd-fitResidSPXLtd$estimate[1])/fitResidSPXLtd$estimate[2]
nSPXL=length(normResidSPXL)
x=qt((1:nSPXL)/(nSPXL+1),df=fitResidSPXLtd$estimate[3])
qqplot(x,sort(normResidSPXL),xlab="t",ylab="residuals", main = "SPXL t QQ Plot Standardized Residuals")
qqline(normResidSPXL)

##Pre-copula transformations and Uhat plot##
uMCD=pt(normResidMCD, df=fitResidMCDtd$estimate[3])
uSPXL=pt(normResidSPXL, df=fitResidSPXLtd$estimate[3])
Uhat=cbind(uMCD,uSPXL)
plot(Uhat, main="MCD Cum. Distribution vs. SPXL Cum. Distribution")

##Fitting copula family to data##
tau=as.numeric(cor.test(uMCD,uSPXL,method="kendall")$estimate)
omega=sin(tau*pi/2)
#t Copula#
Ct = fitCopula(copula=tCopula(dim=2),data=Uhat,method="ml",start=c(omega,6))
llt = loglikCopula(param=Ct@estimate,u=Uhat,copula=tCopula(dim=2))
aict = -2*llt + 2*length(Ct@estimate)
#Clayton Copula#
Cfc = fitCopula(copula=claytonCopula(1,dim=2),data=Uhat,method="ml")
llc = loglikCopula(param=Cfc@estimate,u=Uhat,copula=claytonCopula(dim=2))
aicc = -2*llc + 2*length(Cfc@estimate)
#Frank Copula#
Cfr = fitCopula(copula=frankCopula(1,dim=2),data=Uhat,method="ml")
llf = loglikCopula(param=Cfr@estimate,u=Uhat,copula=frankCopula(dim=2))
aicf = -2*llf + 2*length(Cfr@estimate)
#Gaussian Copula#
Cgauss = fitCopula(copula=normalCopula(dim=2),data=Uhat,method="ml",start=c(omega))
llg = loglikCopula(param=Cgauss@estimate,u=Uhat,copula=normalCopula(dim=2))
aicg = -2*llg + 2*length(Cgauss@estimate)

aict
aicc
aicf
aicg

##Pulling random sample from copula, generate and plot simulated standardized errors from sample## 
X=rCopula(10000, tCopula(dim=2, Ct@estimate[1], df=Ct@estimate[2]))

simulResidMCD=(qt(X[,1], df=fitResidMCDtd$estimate[3])*fitResidMCDtd$estimate[2])+fitResidMCDtd$estimate[1]
simulResidSPXL=(qt(X[,1], df=fitResidSPXLtd$estimate[3])*fitResidSPXLtd$estimate[2])+fitResidSPXLtd$estimate[1]

simulhat=cbind(simulResidMCD, simulResidSPXL)
plot(simulhat, main="MCD & SPXL Simulated Stand. Residual Plot ")

##Finding Sigma for time t+1##
sigmat0MCD = fitMCDtd@sigma.t[length(MCDrtd)]
restMCD=MCDrtd[length(MCDrtd)]
omegaMCD=fitMCDtd@fit$coef[3]
alphaMCD=fitMCDtd@fit$coef[4]
betaMCD=fitMCDtd@fit$coef[5]
sigmat1MCD=sqrt(omegaMCD+(alphaMCD*(restMCD^2))+(betaMCD*(sigmat0MCD^2)))

sigmat0SPXL=fitSPXLtd@sigma.t[length(SPXLrtd)]
restSPXL=SPXLrtd[length(SPXLrtd)]
omegaSPXL=fitSPXLtd@fit$coef[3]
alphaSPXL=fitSPXLtd@fit$coef[4]
betaSPXL=fitSPXLtd@fit$coef[5]
sigmat1SPXL=sqrt(omegaSPXL+(alphaSPXL*(restSPXL^2))+(betaSPXL*(sigmat0SPXL^2)))

##Finding log return for time t+1##
muMCD=fitMCDtd@fit$coef[1]
phiMCD=fitMCDtd@fit$coef[2]
MCDt1=as.numeric(MCD[length(MCD)])
MCDt2=muMCD+(phiMCD*(MCDt1-muMCD))+(sigmat1MCD*simulResidMCD)

muSPXL=fitSPXLtd@fit$coef[1]
phiSPXL=fitSPXLtd@fit$coef[2]
SPXLt1=as.numeric(SPXL[length(SPXL)])
SPXLt2=muSPXL+(phiSPXL*(SPXLt1-muSPXL))+(sigmat1SPXL*simulResidSPXL)

##Summing of simulated log returns##
SumSimulLR=MCDt2+SPXLt2

##Finding VaR, exp shortfall##
VaR=quantile(SumSimulLR, p=.01)
tails=SumSimulLR[(SumSimulLR<=VaR)]
es=mean(tails)
VaR
SumSimulLR if<VaR

##Sums of log returns density plot##
plot(density(SumSimulLR), main="Sum MCD & SPXL Log Return Density")
