diff.var.cycle.poly = function(ts1,ts2,p.order,nharm,npoly,period=12,test.equal.intercept=FALSE,first.step=c(1,1),alpha=0.05,monte.carlo=FALSE,lag.lbtest=10) {
### TESTS EQUALITY OF VAR(P) PLUS ANNUAL/DIURNAL CYCLE MODEL
### HIERARCHY: EQUAL VARIANCE; EQUAL AR; EQUAL CYCLE; EQUAL POLY; EQUAL INTERCEPT
### INPUT:
###		TS1[NTIME1,NSPACE]: TIME SERIES ARRAY 1
###     TS2[NTIME2,NSPACE]: TIME SERIES ARRAY 2 (MUST HAVE SAME SPATIAL DIMENSION)
###     P.ORDER: ORDER OF THE VAR MODEL (can be 0)
###     NHARM: NUMBER OF HARMONICS (can be 0)
###     NPOLY: ORDER OF THE POLYNOMIAL (can be 0)
###		PERIOD: PERIOD OF THE SINUSOIDS (default= 12 for annual cycle with monthly data)
###		TEST.EQUAL.INTERCEPT: SET TO TRUE TO TEST EQUALITY OF INTERCEPTS; OTHERWISE SET TO FALSE
###		FIRST.STEP[2]: THE PHASE OF THE FIRST STEP OF THE ANNUAL/DIURNAL CYCLES
###		ALPHA: SIGNIFICANCE LEVEL
### OUTPUT:

######################################################
########## DEFINE METADATA
######################################################
ts1    = as.matrix(ts1)
ts2    = as.matrix(ts2)
nspace = dim(ts1)[2]

if (nspace != dim(ts2)[2]) stop("ts1 and ts2 should have same number of columns")

alpha.extra = alpha * c(1,1/5)

######################################################
########## SET UP BREAK POINTS AND DEVIANCE LABELS
######################################################
dev.row.names = 'D[0:1]; equal variance'
x.break       = NULL
step          = 0
if (p.order != 0 ) {
	step          = step+1
	x.break       = c(x.break,p.order*nspace)
	dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,']; equal AR model',sep=''))
}
if (nharm != 0 ) {
	step          = step+1
	x.break       = c(x.break,2*nharm)
	dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,']; equal cycle',sep=''))
}
if (npoly != 0 ) {
	step          = step+1
	x.break       = c(x.break,npoly)
	dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,']; equal poly',sep=''))
}
if (test.equal.intercept) {
	step          = step+1
	x.break       = c(x.break,1)
	dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,']; intercept',sep=''))
}
dev.row.names = c(dev.row.names,paste('D[0:',length(dev.row.names),']; total',sep=''))
	
######################################################
########## SET UP TIME SERIES 1
######################################################
ts1.lhs = NULL
ts1.lag = NULL
ts1.cyc = NULL
ts1.pol = NULL
ts1.jvc = NULL
for (ns in 1:nspace) {
	ts1.list  = timeseries2ar.cycle.poly(ts1[,ns],p.order,nharm,npoly,period,first.step=first.step[1])
	if (p.order > 0) ts1.lag = cbind(ts1.lag,as.matrix(ts1.list$y.lag))
	ts1.lhs  = cbind(ts1.lhs,as.matrix(ts1.list$y.lhs))
}
if (nharm   > 0) ts1.cyc = cbind(ts1.cyc,as.matrix(ts1.list$y.cyc))
if (npoly   > 0) ts1.pol = cbind(ts1.pol,as.matrix(ts1.list$y.pol))
ts1.jvc  = cbind(ts1.jvc,as.matrix(ts1.list$jvec))


######################################################
########## SET UP TIME SERIES 2
######################################################
ts2.lhs = NULL
ts2.lag = NULL
ts2.cyc = NULL
ts2.pol = NULL
ts2.jvc = NULL
for (ns in 1:nspace) {
	ts2.list  = timeseries2ar.cycle.poly(ts2[,ns],p.order,nharm,npoly,period,first.step=first.step[2])
	if (p.order > 0) ts2.lag = cbind(ts2.lag,as.matrix(ts2.list$y.lag))
	ts2.lhs  = cbind(ts2.lhs,as.matrix(ts2.list$y.lhs))
}
if (nharm   > 0) ts2.cyc = cbind(ts2.cyc,as.matrix(ts2.list$y.cyc))
if (npoly   > 0) ts2.pol = cbind(ts2.pol,as.matrix(ts2.list$y.pol))
ts2.jvc  = cbind(ts2.jvc,as.matrix(ts2.list$jvec))

######################################################
########## RE-SHAPE LAG DATA
######################################################
if (nspace > 1 & p.order > 0) {
	ntime1       = dim(ts1.lag)[1]
	ntime2       = dim(ts2.lag)[1]
	dim(ts1.lag) = c(ntime1,p.order,nspace)
	dim(ts2.lag) = c(ntime2,p.order,nspace)
	ts1.lag      = aperm(ts1.lag,c(1,3,2))
	ts2.lag      = aperm(ts2.lag,c(1,3,2))
	dim(ts1.lag) = c(ntime1,nspace*p.order)
	dim(ts2.lag) = c(ntime2,nspace*p.order)
}



nest.list = diff.regression.nested.mult(ts1.lhs,ts2.lhs,
   cbind(ts1.lag,ts1.cyc,ts1.pol,ts1.jvc),
   cbind(ts2.lag,ts2.cyc,ts2.pol,ts2.jvc),
   x.break=x.break,alpha=alpha,alpha.extra=alpha.extra,monte.carlo=monte.carlo)

rownames(nest.list$dev.table)       = dev.row.names
rownames(nest.list$dev.table.extra) = dev.row.names

diagnose.list = diagnose.var.cycle.poly(nest.list,x.break,nspace,ts1.cyc,ts1.pol)
nest.list$diagnose.list = diagnose.list

###### PORTMANTAU TEST
lbtest1 = LjungBox(residuals(nest.list$lm1),lags=lag.lbtest,order=p.order)
lbtest2 = LjungBox(residuals(nest.list$lm2),lags=lag.lbtest,order=p.order)
nest.list$lbtest1 = lbtest1
nest.list$lbtest2 = lbtest2


############################################
### COMPUTE ANNUAL CYCLE FORCING (NOT RESPONSE)
### FIT SEPARATELY ON TS1 AND TS2, OTHERWISE "ANNUAL CYCLE" DEPENDS ON BOTH MODELS
### ASSUME COS = 1 IS "JANUARY"
############################################
if (nharm > 0) {
	npic.cyc1    = p.order * nspace + 1:(2*nharm)
	icyc         = min(which(ts1.cyc[,'cos1'] == 1)) + 1:12 - 1
	if (rownames(nest.list$lm1$coef[npic.cyc1,])[1] != 'xmatcos1') stop('accessing wrong beta coefficient?')
	ann.cyc1     = ts1.cyc[icyc,] %*% nest.list$lm1$coef[npic.cyc1,]
	ann.cyc2     = ts1.cyc[icyc,] %*% nest.list$lm2$coef[npic.cyc1,]
	
	nest.list$ann.cyc1 = ann.cyc1
	nest.list$ann.cyc2 = ann.cyc2
} else {
	nest.list$ann.cyc1 = NA
	nest.list$ann.cyc2 = NA
}


nest.list


}




