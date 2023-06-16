timeseries2ar.cycle.poly = function(y,order.pic,nharm,npoly,period=12,first.step=1) {
##### GIVEN TIMESERIES Y, CREATE MATRICES FOR AR(p) PLUS ANNUAL CYCLE MODEL
##### INPUT:
##### 	Y[NTOT]:	TIME SERIES OF LENGTH NTOT
#####   ORDER.PIC:	ORDER OF THE AUTOREGRESSIVE MODEL
#####   NHARM:		NUMBER OF HARMONICS
#####	NPOLY:		ORDER OF POLYNOMIAL CURVE
#####	PERIOD:		PERIOD OF THE CYCLE (DEFAULT = 12 FOR ANNUAL CYCLE IN MONTHLY DATA)
##### 	FIRST.STEP: INTEGER INDICATING THE PHASE OF THE FIRST TIME STEP
##### OUTPUT:LIST$
#####	Y.LHS[NTOT-ORDER.PIC]: TIME SERIES AFTER OMITTING FIRST ORDER.PIC STEPS
#####	Y.LAG[NTOT-ORDER.PIC,ORDER.PIC]: LAGGED PREDICTOR MATRIX FOR Y; EQUALS NULL IF ORDER.PIC = 0
#####	Y.CYC[NTOT-ORDER.PIC,2*NHARM]:    COSINE/SINE PREDICTOR MATRIX; EQUALS NULL IF NHARM = 0
#####	JVEC [NTOT-ORDER.PIC]: VECTOR OF ONES (FOR THE INTERCEPT)

ntot = length(y)
npic = (1+order.pic):ntot
ncyc = (npic - 1) + (first.step - 1)

y.lhs = y[npic]

y.lag = NULL
if (order.pic > 0) {
	for (lag in 1:order.pic) y.lag = cbind(y.lag,y[npic-lag])
	colnames(y.lag) = paste('lag',1:order.pic,sep='')
}

y.cyc = NULL
if (nharm > 0) {
	for (nh in 1:nharm) y.cyc = cbind(y.cyc,cos(2*pi*ncyc*nh/period),sin(2*pi*ncyc*nh/period))
	colnames(y.cyc) = paste(c('cos','sin'),rep(1:nharm,each=2),sep='')
}

y.pol = NULL
if (npoly > 0) {
	y.pol = poly(npic,npoly,simple=TRUE)
	colnames(y.pol) = paste('poly',1:npoly,sep='')
}


jvec = rep(1,length(y.lhs))

list(y.lhs = y.lhs, y.lag = y.lag, y.cyc = y.cyc, y.pol = y.pol, jvec = jvec)
	
}