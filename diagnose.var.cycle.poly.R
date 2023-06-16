diagnose.var.cycle.poly = function(nest.list,x.break,nspace,ts1.cyc,ts1.pol) {

############################################
### DIAGNOSE DISSIMILARITY OF X.MATRIX %*% REGRESSION COEFFICIENTS
############################################
ktot              = nest.list$ktot
ktot.star         = nest.list$ktot.star
ncomm             = nest.list$ncomm
mmat              = array(NA,dim=c(ktot+ktot.star,ktot+ktot.star,ncomm))

### SET UP A-MATRIX AND C-MATRIX
for (nc in 2:ncomm) {
	amat          = NULL
	num.rows      = x.break[nc-1]
	
	if (nc != 2) amat = array(0,dim=c(num.rows,sum(x.break[2:(nc-1)-1])))
	amat          = cbind(amat,diag(1,nrow=num.rows,ncol=x.break[nc-1]))
	
	amat.end1     = NULL
	amat.end2     = NULL
	if (ktot      - dim(amat)[2] > 0) amat.end1 = array(0,dim=c(num.rows,ktot      - dim(amat)[2]))
	if (ktot.star - dim(amat)[2] > 0) amat.end2 = array(0,dim=c(num.rows,ktot.star - dim(amat)[2]))
	
	amat          = cbind(amat,amat.end1,-amat,amat.end2)
	if (nc == 2) amat.list = list(NA) 
	if (nc == 2) cmat.list = list(NA); cmat.list[2] = list(NA)
	amat.list[nc] = list(amat)
	
	if (nc > 2) cmat.list[nc] = list(amat)
	if (nc > 2) amat.list[nc] = list(rbind(amat.list[[nc-1]],amat))
}

xtx.omega1.inv    = chol2inv(chol(t(nest.list$xmat.omega1) %*% nest.list$xmat.omega1))

diagnose.list     = list(list(NA))

for (nc in 2:ncomm) {
	gamma.i           = nest.list$qmat.omega[,,nc-1]/ nest.list$dof.sum
	big.delta         = (nest.list$qmat.omega[,,nc] - nest.list$qmat.omega[,,nc-1]) / nest.list$dof.sum
	
	gamma.eigen       = eigen(gamma.i)
	gamma.sqrt        = gamma.eigen$vectors %*% diag(  sqrt(gamma.eigen$values),nrow=nspace,ncol=nspace) %*% t(gamma.eigen$vectors)
	gamma.inv.sqrt    = gamma.eigen$vectors %*% diag(1/sqrt(gamma.eigen$values),nrow=nspace,ncol=nspace) %*% t(gamma.eigen$vectors)

	mmat[,,nc] = diag(1,nrow=ktot+ktot.star,ncol=ktot+ktot.star) - xtx.omega1.inv %*% 
	   t(amat.list[[nc]]) %*% chol2inv(chol(amat.list[[nc]] %*% xtx.omega1.inv %*% t(amat.list[[nc]]))) %*% amat.list[[nc]] 
	
	if (nc == 2) {
		lit.delta    = amat.list[[nc]] %*% nest.list$lm.omega[[1]]$coef
		theta        = chol2inv(chol(amat.list[[nc]] %*% xtx.omega1.inv %*% t(amat.list[[nc]]))) / nest.list$dof.sum
		delta.check  = t(lit.delta) %*% theta %*% lit.delta
		# if (!isTRUE(all.equal(delta.check,big.delta))) stop('deltas inconsistent')	

		beta.stack = mmat[,,nc] %*% nest.list$lm.omega[[1]]$coef
		beta1      = nest.list$lm.omega[[1]]$coef[1:x.break[1],]
		beta2      = nest.list$lm.omega[[1]]$coef[1:x.break[1]+ktot,]
		# if (!isTRUE(all.equal(as.numeric(beta.stack[1:ktot,]),as.numeric(nest.list$lm.omega[[nc]]$coef[1:ktot,])))) stop('betas inconsistent')
		# if (!isTRUE(all.equal(as.numeric(beta1-beta2),as.numeric(lit.delta)))) stop('beta1/beta2 inconsistent with lit.delta')
		
	} else {
		gmat       = amat.list[[nc-1]] %*% xtx.omega1.inv %*% t(amat.list[[nc-1]])
		hmat       = amat.list[[nc-1]] %*% xtx.omega1.inv %*% t(cmat.list[[nc  ]])
		jmat       = cmat.list[[nc  ]] %*% xtx.omega1.inv %*% t(cmat.list[[nc  ]])
		lmat       = chol2inv(chol(jmat - t(hmat) %*% chol2inv(chol(gmat)) %*% hmat))
		zmat       = xtx.omega1.inv %*% t(mmat[,,nc-1]) %*% t(cmat.list[[nc]]) %*% lmat
		mmat.check = mmat[,,nc-1] - zmat %*% cmat.list[[nc  ]] %*% mmat[,,nc-1]
		# if (!isTRUE(all.equal(mmat[,,nc],mmat.check))) stop('m-matrix inconsistent')
		
		beta.stack  = mmat[,,nc-1] %*% nest.list$lm.omega[[1]]$coef
		
		lit.delta   = cmat.list[[nc  ]] %*% beta.stack
		theta       = lmat / nest.list$dof.sum
		# delta.check = t(lit.delta) %*% t(zmat) %*% t(nest.list$xmat.omega1) %*% nest.list$xmat.omega1 %*% zmat %*% lit.delta / nest.list$dof.sum
		# if (!isTRUE(all.equal(big.delta,delta.check))) stop('delta 1 inconsistent')
		# if (!isTRUE(all.equal(big.delta,t(lit.delta) %*% theta %*% lit.delta))) stop('delta 2 inconsistent')
		# theta.check = t(zmat) %*% t(nest.list$xmat.omega1) %*% nest.list$xmat.omega1 %*% zmat / nest.list$dof.sum
		# if (!isTRUE(all.equal(theta.check,theta))) stop('theta inconsistent')
		
		kpic        = sum(x.break[1:(nc-2)]) + 1:x.break[nc-1]
		beta1       = beta.stack[kpic,]
		beta2       = beta.stack[kpic+ktot,]
		# if (!isTRUE(all.equal(as.numeric(beta1-beta2),as.numeric(lit.delta)))) stop('beta1/beta2 inconsistent with lit.delta')
	}	
	
	theta.eigen    = eigen(theta)
	theta.sqrt     = theta.eigen$vectors %*% diag(  sqrt(theta.eigen$values),nrow=length(theta.eigen$values),ncol=length(theta.eigen$values)) %*% t(theta.eigen$vectors)
	theta.inv.sqrt = theta.eigen$vectors %*% diag(1/sqrt(theta.eigen$values),nrow=length(theta.eigen$values),ncol=length(theta.eigen$values)) %*% t(theta.eigen$vectors)
	phi.mat        = theta.sqrt %*% lit.delta %*% gamma.inv.sqrt
	phi.svd        = svd(phi.mat)
	qx.mat         = theta.inv.sqrt %*% phi.svd$u
	px.mat         = theta.sqrt     %*% phi.svd$u
	qy.mat         = gamma.inv.sqrt %*% phi.svd$v
	py.mat         = gamma.sqrt     %*% phi.svd$v
	
	py1.mat        = t(px.mat) %*% beta1
	py2.mat        = t(px.mat) %*% beta2

	# if (!isTRUE(all.equal(py1.mat-py2.mat,phi.svd$d * t(py.mat)))) stop('beta1/beta2 inconstent propagations')	
	# if (!isTRUE(all.equal(t(px.mat) %*% lit.delta, phi.svd$d * t(py.mat)))) stop('pxT delta inconsistent')
	# if (!isTRUE(all.equal(big.delta,py.mat %*% diag(phi.svd$d^2) %*% t(py.mat)))) stop('delta 3 inconsistent')
	# if (!isTRUE(all.equal(phi.svd$d^2,nest.list$cda.svalues[,nc]))) stop('s. values inconsistent')
		
	if (nc == 2) {
		qx.dot.mat  = NA
		px.dot.mat  = NA
		qx1.dot.mat = NA
		qx2.dot.mat = NA
	} else {
		if (nc == 3) {
			### ANNUAL CYCLE SUBSET
			icyc        = min(which(ts1.cyc[,'cos1'] == 1)) + 1:12 - 1
			x.sub       = ts1.cyc[icyc,]
		} else {
			print('have not defined polynomial X.SUB')
			x.sub       = nest.list$xmat.omega1[,x.break[nc-1-1]+1:x.break[nc-1]]
		}
		
		xtx.sub.inv = chol2inv(chol(t(x.sub) %*% x.sub))
		qx.dot.mat  = x.sub %*% qx.mat
		px.dot.mat  = x.sub %*% xtx.sub.inv %*% px.mat
	
		kpic1           = sum(x.break[2:(nc-1)-1])+1:x.break[nc-1]
		kpic2           = kpic1+ktot
		x.b1            = x.sub %*% beta.stack[kpic1,]
		x.b2            = x.sub %*% beta.stack[kpic2,]
		# if (!isTRUE(all.equal(x.b1-x.b2,x.sub %*% lit.delta))) stop('difference is not equal to difference in betas')
		
		qx1.dot.mat     = x.b1 %*% qy.mat
		qx2.dot.mat     = x.b2 %*% qy.mat
		# if (!isTRUE(all.equal(qx1.dot.mat-qx2.dot.mat,qx.dot.mat %*% diag(phi.svd$d)))) stop("qxdot inconsistent")	
	}

	diagnose.list[nc] = list(list(s.values = phi.svd$d, dev=nest.list$dof.sum *log(1+phi.svd$d^2),
                         qx.mat = qx.mat, qy.mat = qy.mat, px.mat = px.mat, py.mat = py.mat,
                         qx1.dot.mat = qx1.dot.mat, qx2.dot.mat = qx2.dot.mat, 
                         qx.dot.mat = qx.dot.mat, px.dot.mat = px.dot.mat,
                         py1.mat = py1.mat, py2.mat = py2.mat))
	
	

	######  PLOTS
	# xrange = c(1,12)
	# yrange = range(qx1.dot.mat[,1],qx2.dot.mat[,1])
	# par(mfcol=c(2,1),mar=c(5,5,3,1))
	# plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='month',ylab='')
	# lines(qx1.dot.mat[1:xrange[2],1],lwd=2)
	# lines(qx2.dot.mat[1:xrange[2],1],lwd=2,col='red')
	
	# py.space = laplacian.list$eof[,1:nspace] %*% py.mat[,1]
	# plot_latlon_v4(laplacian.list$lon,laplacian.list$lat,py.space,shrinkdomain=TRUE)



}

names(diagnose.list) = paste('Omega',1:ncomm,sep='')

diagnose.list
	
}