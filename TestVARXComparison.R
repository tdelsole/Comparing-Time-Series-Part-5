rm(list=ls())

### APPLY VARX(P) COMPARISON CODES


nsamp     = 100		## sample size
nspace    = 4		## dimension of VARX model
phi.base  = 0.8		## lag-1 correlation of most persistent variable
nspin     = 100		## length of spinup for VARX model
nharm     = 2		## number of harmonics of periodic cycle
npoly     = 0		## order of the polynomial
p.order   = 1		## p: order of VARX model (only p=1 is allowed here)
constant  = 0		## intercept term in VARX model
n.period  = 12		## period of the fundamental harmonic

set.seed(4)

# dir.portes  = '/Users/delsole/R/portes/R/'
# dir.Rlib    = '/Users/delsole/R/delsole_tools/'
dir.portes  = NULL
dir.Rlib    = NULL
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.portes,'LjungBox.R',sep=''))

#### MULTIVARIATE
source(paste(dir.Rlib,'diff.var.cycle.poly.R',sep=''))
source(paste(dir.Rlib,'diff.regression.nested.mult.R',sep=''))
source(paste(dir.Rlib,'timeseries2ar.cycle.poly.R',sep=''))
source(paste(dir.Rlib,'diagnose.var.cycle.poly.R',sep=''))


#############################################
##########  STOCHASTIC MODEL PARAMETERS
##########  NOISE COVARIANCE MATRIX = IDENTITY; 
##########  PHI_1 = UPPER TRIANGULAR WITH DIAGONAL ELEMENTS < 1
##########  RANDOM ANNUAL CYCLE COEFFICIENTS
#############################################
if (nharm > 0) {
	time.step    = 1:nsamp - 1
	cos.sin.coef = array(rnorm(2*nharm*nspace),dim=c(2*nharm,nspace))
	fmat         = NULL
	for (nh in 1:nharm) fmat = cbind(fmat,cos(2*pi*nh*time.step/n.period),sin(2*pi*nh*time.step/n.period))
} else {
	cos.sin.coef = 0
	fmat = 0
}


phi.mat   = array(0,dim=c(nspace,nspace))
for (j in 1:nspace) for (i in 1:j) if (i == j) phi.mat[i,j] = phi.base^i else phi.mat[i,j] = rnorm(1)

#############################################
##########  GENERATE DATA
#############################################
spin1     = array(NA,dim=c(nspace,nspin))
spin2     = array(NA,dim=c(nspace,nspin))
spin1[,1:p.order] = 0
spin2[,1:p.order] = 0
for (n in (1+p.order):nspin) spin1[,n] = phi.mat %*% spin1[,n-1] + rnorm(nspace) + constant + as.numeric(fmat[n,] %*% cos.sin.coef)
for (n in (1+p.order):nspin) spin2[,n] = phi.mat %*% spin2[,n-1] + rnorm(nspace) + constant + as.numeric(fmat[n,] %*% cos.sin.coef)

ts1       = array(NA,dim=c(nspace,nsamp))    
ts2       = array(NA,dim=c(nspace,nsamp))    
ts1[,1:p.order] = spin1[,nspin+1:p.order-p.order]
ts2[,1:p.order] = spin2[,nspin+1:p.order-p.order]
for (n in (1+p.order):nsamp) ts1[,n] = phi.mat %*% ts1[,n-1] + rnorm(nspace) + constant + as.numeric(fmat[n,] %*% cos.sin.coef)
for (n in (1+p.order):nsamp) ts2[,n] = phi.mat %*% ts2[,n-1] + rnorm(nspace) + constant + as.numeric(fmat[n,] %*% cos.sin.coef)

ts1 = t(ts1)
ts2 = t(ts2)


#############################################
##########  PLOT
#############################################
yrange = range(ts1,ts2)
par(mfrow=c(2,2),mar=c(5,5,3,1))
for (ns in 1:4) {
	plot(ts1[, ns],type='l',ylim=yrange,xlab='time',ylab=ns)
	lines(ts2[, ns],col='red')
}

#############################################
##########  COMPUTE DEVIANCES WITH NEW CODE
#############################################
# stop()
# test.equal.intercept=FALSE; first.step=c(1,1); alpha=0.05; monte.carlo=FALSE; lag.lbtest=10
diff.new = diff.var.cycle.poly(ts1,ts2,p.order,nharm,npoly,monte.carlo=TRUE)

print(diff.new$dev.table)

#############################################
##########  PLOT RESULTS
#############################################
par(mfrow=c(2,2),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
col.pic = c('skyblue2','green2','pink2')

#### BARPLOT FOR SUBDEVIANCES
sub.devs = diff.new$dev.table[1:3,'deviance']
names(sub.devs) = c('noise','AR','cycle')
barplot(sub.devs,las=2,col=col.pic,ylab='sub-deviance')

### DISCRIMINANTS OF NOISE (I.E., WHITENED VARIANCE)
for (np in 1:3) {
	if (np == 1) {
		plot.vals  = diff.new$cda.all[[np]]$evalues
		plot.title = 'Noise Variance Ratio'
		plot.line.col = col.pic[np]
		plot.ylab     = 'Noise Variance Ratio'
		plot.crit     = c(diff.new$cda.noise.ratio.min,diff.new$cda.noise.ratio.max)
	} else if (np == 2) {
		plot.vals  = diff.new$cda.all[[np]]$devs
		plot.title = 'AR Sub-Deviance'
		plot.line.col = col.pic[np]
		plot.ylab     = 'AR Sub-Deviance'
		plot.crit     = diff.new$dev.crit.paramdiff[np]	
	} else if (np == 3) {
		plot.vals  = diff.new$cda.all[[np]]$devs
		plot.title = 'Cycle Sub-Deviance'
		plot.line.col = col.pic[np]
		plot.ylab     = 'cycle sub-deviance'				
		plot.crit     = diff.new$dev.crit.paramdiff[np]	
	}
	yrange = range(plot.vals,plot.crit)
	plot(plot.vals,type='l',lwd=2,col=plot.line.col,xaxt='n',xlab='discriminant',ylab=plot.ylab,ylim=yrange)
	axis(1,at=1:nspace)
	points(plot.vals,pch=19,col=plot.line.col)
	title(main=plot.title,line=0.5)
	abline(h=plot.crit,lty='dashed',col=plot.line.col)
}

