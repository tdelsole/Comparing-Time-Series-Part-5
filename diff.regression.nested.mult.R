diff.regression.nested.mult = function(yvec,yvec.star,xmat,xmat.star,x.break,alpha=0.05,monte.carlo=FALSE,ntrials=1000,alpha.extra=c(0.05,0.01)) {
###### TEST EQUALITY OF REGRESSION MODELS THROUGH NESTED HYPOTHESES:
###### YVEC = X2  B2  + ... + XL  BL  + E
###### YVEC = X2* B2* + ... + XL* BL* + E*
###### IMPORTANT: USER SHOULD INCLUDE INTERCEPT IN XMAT AND XMAT.STAR
###### HYPOTHESES: 
######	OMEGA.0: UNRESTRICTED
######  OMEGA.1: EQUALITY OF NOISE VARIANCE (VAR[E] = VAR[E*]):  X.COMM[1] = 0
######  OMEGA.2: OMEGA.1 AND B2 = B2*; LENGTH(B2) = X.BREAK[1] = X.COMM[2]
######  OMEGA.3: OMEGA.2 AND B3 = B3*; LENGTH(B3) = X.BREAK[2] = X.COMM[3]
######  ...
###### 	OMEGA.L: OMEGA.L-1 AND BL = BL*; TESTED ONLY IF X.BREAK[L-1] IS SPECIFIED
###### INPUT:
###### 	YVEC     [NTOT,NSPACE]: 	Y      TIMESERIES
###### 	YVEC.STAR[NTOT.STAR,NSPACE]:Y.STAR TIMESERIES
######	XMAT     [NTOT,KTOT]:		MATRIX [X2 ,X3 ,...,XL ]
######	XMAT.STAR[NTOT,KTOT.STAR]:	MATRIX [X2*,X3*,...,XL*]
######  X.BREAK  [LTOT]:			NUMBER OF COLUMNS IN X2, X3, ... XL; SUM(X.BREAK) <= KTOT
######	ALPHA:						SIGNIFICANCE LEVEL; DEFAULT IS 5%
######  MONTE.CARLO[LOGICAL]:		SET TO 'TRUE' TO DO MONTE CARLO ESTIMATES OF CRITICAL VALUES (EXPENSIVE!)
######	NTRIALS:					NUMBER OF MONTE CARLO TRIALS FOR ESTIMATING DEVIANCE [OMEGA.0 - OMEGA.1]
######  ALPHA.EXTRA[NALPHA]:		ADDITIONAL SIGNIFICANCE LEVELS TO ASSESS; DEFAULT IS 5% AND 1%
###### OUTPUT:
######	 DEV.TABLE[LTOT+1,7]: ROWS: LTOT SUB-DEVIANCES, PLUS TOTAL DEVIANCE.  COLUMNS: DEV, CRITICAL (ASYM), PVAL (ASYM), ALPHA, CRITICAL (MC)
######	 LM.OMEGA[[LTOT]]: OUTPUT CAPTURED FROM 'LM' FOR HYPOTHESES OMEGA_1,...,OMEGA_LTOT
###### 	LM1: OUTPUT CAPTURED FROM 'LM' FOR FITTING MODEL 1 ALONE
###### 	LM2: OUTPUT CAPTURED FROM 'LM' FOR FITTING MODEL 2 ALONE
###### 	ALPHA.STEP: THE INDIVIDUAL ALPHA LEVEL TO CONTROL FAMILYWISE ERROR RATE
###### 	DEV.TABLE.EXTRA[LTOT+1,4]: SIGNIFICANCE LEVELS FOR LTOT SUB-DEVIANCES, PLUS TOTAL DEVIANCE, FOR ALPHA AND ALPHA/5, ASYMPTOTIC AND MC.   
###### 	QMAT1, QMAT2, QMAT.OMEGA, XMAT.OMEGA1: ARRAYS USED IN DIAGNOSIS
###### 	CDA.NOISE.RATIO.MAX: ALPHA/2 SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT NOISE VARIANCE RATIO 
###### 	CDA.NOISE.RATIO.MIN: ALPHA/2 SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT NOISE VARIANCE RATIO 
###### 	DEV.CRIT.PARAMDIFF[LTOT]: SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT FOR D[I:I+1].  I=1 IS NA BECAUSE ONLY RATIO IS REPORTED
###### 	CDA.ALL[[LTOT]]: LIST OF (EIGENVALUES, DEVIANCES, Q, P) FROM CDA OF D[I:I+1]
###### 	KTOT: TOTAL NUMBER OF PREDICTORS IN FIRST MODEL
###### 	KTOT.STAR: TOTAL NUMBER OF PREDICTORS IN 2ND MODEL
###### 	NCOMM: NUMBER OF PREDICTOR BLOCKS TO BE TESTED FOR EQUALITY ("KAPPA")
###### 	DOF1: DEGREES OF FREEDOM FOR FITTING 1ST MODEL
###### 	DOF2: DEGREES OF FREEDOM FOR FITTING 2ND MODEL
###### 	DOF.SUM: DOF1 + DOF2


##############################
###### METADATA
##############################
x.comm    = c(0,x.break) # append 0 at the start, corresponding to no common predictors
xmat      = as.matrix(xmat)
xmat.star = as.matrix(xmat.star)
yvec      = as.matrix(yvec)
yvec.star = as.matrix(yvec.star)

ncomm     = length(x.comm)
ntot      = dim(yvec)[1]
ntot.star = dim(yvec.star)[1]
nspace    = dim(yvec)[2]
if (sum(x.comm) > dim(xmat)     [2]) stop('x.comm inconsistent with xmat')
if (sum(x.comm) > dim(xmat.star)[2]) stop('x.comm inconsistent with xmat.star')

alpha.step  = 1 - (1-alpha)^(1/ncomm)
param.num0  = NA

nalpha           = length(alpha.extra)
alpha.extra.step = 1 - (1-alpha.extra)^(1/ncomm)

########################################
####### EQUALITY OF VARIANCES (OMEGA_1)
########################################
lm1     = lm(yvec      ~ xmat       - 1)
lm2     = lm(yvec.star ~ xmat.star  - 1)
qmat1   = t(residuals(lm1)) %*% residuals(lm1) 
qmat2   = t(residuals(lm2)) %*% residuals(lm2) 
dof1    = lm1$df.residual
dof2    = lm2$df.residual
dof.sum = dof1 + dof2

########################################
####### COMPUTE COVARIANCE MATRICES FOR NESTED HYPOTHESES
########################################
y.all         = rbind(yvec,yvec.star)
qmat.omega    = array(NA,dim=c(nspace,nspace,ncomm))
dof.omega     = as.numeric(rep(NA,ncomm))
dev.omega     = as.numeric(rep(NA,ncomm))
param.num     = as.numeric(rep(NA,ncomm))
phi           = as.numeric(rep(NA,ncomm))
pval.asym     = as.numeric(rep(NA,ncomm))

for (nb in 1:ncomm) {
	ncol.comm      = sum(x.comm[1:nb])
	ncol.diff      = ncol(xmat)      - ncol.comm
	ncol.diff.star = ncol(xmat.star) - ncol.comm
	
	if (ncol.comm == 0) {
		x.all = NULL
	} else {
		x.all  = rbind(xmat[,1:ncol.comm,drop=FALSE],xmat.star[,1:ncol.comm,drop=FALSE])
	}
	
	if (ncol.diff > 0) {
		x.diff      = cbind(xmat     [,1:ncol.diff      + ncol.comm])
		x.all       = cbind(x.all,rbind(x.diff,array(0,dim=c(ntot.star,ncol.diff.star))))
	}
	if (ncol.diff.star > 0) {
		x.diff.star = cbind(xmat.star[,1:ncol.diff.star + ncol.comm])		
		x.all       = cbind(x.all,rbind(array(0,dim=c(ntot,ncol.diff)),x.diff.star))
	}
	
	lm.nest          = lm(y.all ~ x.all - 1)
	qmat.omega[,,nb] = t(residuals(lm.nest)) %*% residuals(lm.nest)  
	dof.omega[nb]    = lm.nest$df.residual
	if (nb == 1) lm.omega = list(lm.nest) else lm.omega[nb] = list(lm.nest)
	param.num[nb]    = nspace * ncol(x.all) + nspace*(nspace+1)/2
	phi[nb]          = ncol(x.all)
	
	if (ncol.comm == 0) param.num0 = nspace * ncol(x.all) + nspace*(nspace+1)
	if (ncol.comm == 0) xmat.omega1 = x.all
}



#### DEVIANCE 0:1
q.gev           = gev(qmat1/dof1,qmat2/dof2)
dev.omega.check = dof.sum * log(det((dof1*qmat1+dof2*qmat2)/dof.sum)) - dof1 * log(det(qmat1)) - dof2 * log(det(qmat2))
dev.omega[1]    = dof1 * sum(log((dof1+dof2/q.gev$lambda)/dof.sum)) + dof2 * sum(log((dof1*q.gev$lambda + dof2)/dof.sum))
if (abs(dev.omega.check-dev.omega[1]) > 1.e-10) stop('inaccurate noise deviance detected')
devs            = dof1 * log((dof1+dof2/q.gev$lambda)/dof.sum) + dof2 * log((dof1*q.gev$lambda + dof2)/dof.sum)
cda.all         = list(list(evalues = q.gev$lambda, devs = devs, q = q.gev$q, p = qmat2 %*% q.gev$q /dof2))

#### OTHER DEVIANCES
for (nb in 2:ncomm) {
	q.gev         = gev(qmat.omega[,,nb]/dof.sum,qmat.omega[,,nb-1]/dof.sum)
	# dev.check     = dof.sum * log(det(qmat.omega[,,nb])/det(qmat.omega[,,nb-1]))
	dev.omega[nb] = dof.sum * sum(log(q.gev$lambda))
	devs             = dof.sum * log(q.gev$lambda)
	cda.all[nb] = list(list(evalues = q.gev$lambda-1, devs = devs, q = q.gev$q, p = qmat.omega[,,nb-1] %*% q.gev$q /dof.sum))
}

ktot              = dim(xmat     )[2]
ktot.star         = dim(xmat.star)[2]


###################################################
######### MONTE CARLO 
###################################################
dev.crit.mc             = as.numeric(rep(NA,ncomm))
dev.crit.asym           = as.numeric(rep(NA,ncomm))
dev.crit.extra.mc       = array(NA,dim=c(nalpha,ncomm))
dev.crit.extra.asym     = array(NA,dim=c(nalpha,ncomm))
dev.total.crit.mc       = NA
dev.total.crit.extra.mc = NA
cda.noise.ratio.max     = NA
cda.noise.ratio.min     = NA
dev.crit.paramdiff      = array(NA,dim=c(ncomm))
names(dev.crit.paramdiff)  = paste('Omega',1:ncomm,sep='')
names(cda.all)             = paste('Omega',1:ncomm,sep='')


if (monte.carlo) {
	### DEVIANCE FOR OMEGA.0 VS OMEGA.1
	identity      = diag(1,nrow=nspace,ncol=nspace)
	dev.mc        = array(NA,dim=c(ntrials,ncomm))
	cda.noise.ratio.max.mc = array(NA,dim=c(ntrials))
	cda.noise.ratio.min.mc = array(NA,dim=c(ntrials))
	dev.decomp.paramdiff.mc = array(NA,dim=c(ntrials,ncomm))
	for (nt in 1:ntrials) {
		q1.mc = rWishart(1,dof1,identity)[,,1]
		q2.mc = rWishart(1,dof2,identity)[,,1]
		# lambda = gev(q1.mc,q2.mc)$lambda
		# dev.check = sum( dof.sum * log((lambda+1)/dof.sum) - dof1 * log(lambda/dof1) - dof2 * log(1/dof2))
		# dev.check  = dof.sum * log(det((q1.mc + q2.mc)/dof.sum)) - dof1 * log(det(q1.mc/dof1)) - dof2 * log(det(q2.mc/dof2))
		# print(c(dev.mc[nt],dev.check))
		lambda = gev(q1.mc/dof1,q2.mc/dof2)$lambda
		dev.mc[nt,1]  = sum( dof.sum * log((dof1 * lambda + dof2)/dof.sum) - dof1 * log(lambda) )
		cda.noise.ratio.max.mc[nt] = lambda[1]
		cda.noise.ratio.min.mc[nt] = lambda[nspace]
		# print(c(dev.mc[nt,1],dev.check))
		
		q3.old    = q1.mc + q2.mc
		for (nb in 2:ncomm) {
			q3.mc  = rWishart(1,phi[nb-1]-phi[nb],identity)[,,1]
			lambda = gev(q3.mc,q3.old)$lambda
			dev.mc[nt,nb] = dof.sum*sum(log(1+lambda))
			dev.decomp.paramdiff.mc[nt,nb] = dof.sum * log(1+lambda[1])
			# dev.check     = dof.sum*log(det(q1.mc + q2.mc + q3.mc)/det(q1.mc + q2.mc))
			# print(c(dev.mc[nt,nb],dev.check))
			q3.old = q3.old + q3.mc
		}
		
	}
	
	for (nb in 1:ncomm) dev.crit.mc      [ nb] = quantile(dev.mc[,nb],probs=1-alpha.step)
	for (nb in 1:ncomm) dev.crit.extra.mc[,nb] = quantile(dev.mc[,nb],probs=1-alpha.extra.step)
	dev.total.crit.mc       = quantile(rowSums(dev.mc),probs=1-alpha)
	dev.total.crit.extra.mc = quantile(rowSums(dev.mc),probs=1-alpha.extra)
	
	cda.noise.ratio.max = quantile(cda.noise.ratio.max.mc,probs=1-alpha/2)
	cda.noise.ratio.min = quantile(cda.noise.ratio.min.mc,probs=alpha/2)
	for (nb in 2:ncomm) dev.crit.paramdiff[nb] = quantile(dev.decomp.paramdiff.mc[,nb],probs=1-alpha)
	
	
}




# for (nb in 2:ncomm) dev.crit[nb] = 
  # (ntot + ntot.star) * log( 1 + qf(alpha.step,dof.omega[nb]-dof.omega[nb-1],dof.omega[nb-1],lower.tail=FALSE) * (dof.omega[nb]-dof.omega[nb-1])/dof.omega[nb-1])
  
#### ASYMPTOTICS
dev.crit.asym      [ 1] = qchisq(alpha.step      ,param.num0 - param.num[1],lower.tail=FALSE)
dev.crit.extra.asym[,1] = qchisq(alpha.extra.step,param.num0 - param.num[1],lower.tail=FALSE)
for (nb in 2:ncomm) dev.crit.asym      [ nb] = qchisq(alpha.step      ,param.num[nb-1] - param.num[nb],lower.tail=FALSE)
for (nb in 2:ncomm) dev.crit.extra.asym[,nb] = qchisq(alpha.extra.step,param.num[nb-1] - param.num[nb],lower.tail=FALSE)

pval.asym[1] = pchisq(dev.omega[1],param.num0 - param.num[1],lower.tail=FALSE)
for (nb in 2:ncomm) pval.asym[nb] = pchisq(dev.omega[nb],param.num[nb-1] - param.num[nb],lower.tail=FALSE)

#### TOTAL DEVIANCE
dev.total       = sum(dev.omega)
dev.total.check = dof.sum * log(det(as.matrix(qmat.omega[,,ncomm])/dof.sum)) -
    dof1*log(det(as.matrix(qmat1)/dof1)) - dof2*log(det(as.matrix(qmat2)/dof2))
dev.total.pval.asym       = pchisq(dev.total  ,param.num0 - param.num[ncomm],lower.tail=FALSE)
dev.total.crit.asym       = qchisq(alpha      ,param.num0 - param.num[ncomm],lower.tail=FALSE)
dev.total.crit.extra.asym = qchisq(alpha.extra,param.num0 - param.num[ncomm],lower.tail=FALSE)


########################################
#### SUMMARY TABLE
########################################
dev.table = cbind(dev.omega,dev.crit.asym,pval.asym,alpha.step,dev.crit.mc)
dev.table = rbind(dev.table,c(dev.total,dev.total.crit.asym,dev.total.pval.asym,alpha,dev.total.crit.mc))
rownames(dev.table) = c(paste('D[',1:ncomm-1,':',1:ncomm,']',sep=''),paste('D[0:',ncomm,']',sep=''))
colnames(dev.table) = c('deviance','crit(chisq)','pval(chisq)','alpha','crit(MC)')

dev.table.extra = cbind(t(cbind(dev.crit.extra.asym,dev.total.crit.extra.asym)),t(cbind(dev.crit.extra.mc,dev.total.crit.extra.mc)))
rownames(dev.table.extra) = rownames(dev.table)
colnames(dev.table.extra) = c(paste('asym',(1-alpha.extra)*100,'%',sep=''),paste('MC',(1-alpha.extra)*100,'%',sep=''))

########################################
#### OUTPUT THE RESULTS
########################################
list(dev.table = dev.table, lm.omega = lm.omega, lm1=lm1, lm2=lm2, alpha.step=alpha.step, dev.table.extra=dev.table.extra,
  qmat1 = qmat1, qmat2 = qmat2, qmat.omega = qmat.omega, 
  cda.noise.ratio.max = cda.noise.ratio.max, cda.noise.ratio.min = cda.noise.ratio.min, 
  dev.crit.paramdiff = dev.crit.paramdiff, cda.all = cda.all,
  xmat.omega1 = xmat.omega1, ktot = ktot, ktot.star = ktot.star, ncomm = ncomm, dof1 = dof1, dof2 = dof2, dof.sum = dof.sum)

}