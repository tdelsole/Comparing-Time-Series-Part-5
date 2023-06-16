gev <- function (A,B){
# Computes generalized eigenvalue problem A*x=lambda*B*x 
#INPUT: 
# Two matrcies A(ntrun,ntrun) and B(ntrun,ntrun) 

#OUTPUT:
#  q(ntrun,ntrun): projection vectors for the APT-optimals
#  lambda(ntrun) : eigen-values 

# Check if A and B are the same dimension square matrices!
if (length(A) == 1 ) {
	if ( length(B) != 1 ) stop('A and B not dimensioned correctly')
	lambda = A/B
	q      = 1
} else {
	if (!all(dim(A)==dim(B))) stop("dim(A)!=dim(B)")
	if (dim(A)[1]!=dim(A)[2]) stop("A is not square")

    lprint = TRUE
	ntrun = dim(A)[1]
	
	# cholesky decomposition B = U^T U.  
	B.chol = chol(B)
	
	# whitened A = (U^T)^(-1) A U^(-1)
	A.white = forwardsolve(t(B.chol),t(forwardsolve(t(B.chol),A)))
	
	# eigen vectors of whitened A
	eigen.chol = eigen(A.white,symmetric=TRUE)
	
	# Transform eigenvectors evec = U^(-1) evec.white
	q = backsolve(B.chol,eigen.chol$vectors)

	# normalize vectors to qTBq = 1
	amp = vector()
	for ( n in 1:ntrun ) amp[n] = q[,n] %*% B %*% q[,n]
	amp = sqrt(amp)
	for (n in 1:ntrun ) q[,n] = q[,n]/amp[n]
    
    # eigen-values 
    lambda = eigen.chol$values
	
}
  
#    if (lprint) print("Generalized e-value problem finished")
list(lambda=lambda, q=q)
}