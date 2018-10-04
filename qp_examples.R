#
# Example quadratic program data
#

# A QP structured in the style of quadprog
QP <- setClass(
    "QP",
    slots = c("Dmat", "dvec", "Amat", "bvec", "meq"),
    prototype=list(
        meq=0
    )
)

#########################################################
# Random QP Generation
#########################################################
randomQP <- function(N, d){
    # Generate matrices for a random QP
    # 
    # This method returns Dmat an SPD matrix and dvec a 
    # postive random vector defining the quadratic form
    #   -d^T x + 1/2 x^T D x
    # We also return Amat and bvec which give the constraint
    # relation
    #    Ax >= b
    # bvec is a ranod positive random vector of length N/d and Amat is
    # a binary matrix of dimension (N/d) x N.  Amat is defined
    # so that (Ax)_i is the sum of the i through (i+d) components of x.
    #
    # NOTE: differs with quadprog input format by using A rather than A^T
    # in the constraint.
    if(N%%d) stop("N must be divisible by d")
    
    # build an SPD matrix
    G <- matrix(rnorm(N^2, 0, 1), N, N)
    M <- lower.tri(matrix(1, N, N), diag=FALSE)
    G[!M] <- 0
    G <- t(G)%*%G 
    G <- G/sqrt(sum(G^2)) + diag(1,N,N)
    
    # QP data
    qp <- QP(Dmat=G, 
             dvec=1+matrix(runif(N), N, 1), 
             Amat=rbind(kronecker(diag(N/d), matrix(1, 1, d))),
             bvec=rbind(matrix(runif(N), N/d, 1)))
    return(qp)
}

######################################################
# Circus tent example
# 
# The original problem is from the MathWorks MATLAB demo at:
#  http://www.mathworks.com/help/optim/examples/large-scale-bound-constrained-quadratic-programming.html
# 
######################################################
# make a tent pole ensamble; recommeneded to use power of 2 for N to maintain symmetry
makeTentEnsemble <- function(nBase=2**5, nFold=1){
    n <- floor(nBase/2)       # number of rows and number of columns per quandrant
    m <- floor(n/2)       # midpoint row and column index in the first quadrant
    q2 = matrix(0,n,n)
    q2[m:(m+1), m:(m+1)] = .3
    q2[(n-1):n, (n-1):n] = 1/2
    q1 <- q2[,n:1]
    top  <- cbind(q2,q1)
    z <- rbind(top,top[n:1,])
    
    # use nFold to stack copies of z in an nFold by nFold grid
    # rep and stack columns then rep and stack the result down the row
    if(nFold>1) z<-apply(apply(z,2, rep, nFold), 1, rep, nFold)  
    
    # xy axis builder
    dd <- dim(z)
    x <- (1:dd[1])/dd[1]
    y <- (1:dd[2])/dd[2]
    
    # pack result in a list
    res = list(x=x, y=y, z=z)
    return(res)
}

# 2D discrete laplacian
Laplacian2D <- function(lx,ly, nx,ny){
    hx <- lx/(nx-1)
    hy <- ly/(ny-1)
    tr_x <- c(2,-1,rep(0,nx-2))
    tr_y <- c(2,-1,rep(0,ny-2))
    Tx <- toeplitz(tr_x)/(hx^2)
    Ty <- toeplitz(tr_y)/(hy^2)
    Ix <- diag(nx)
    Iy <- diag(ny)
    L <- kronecker(Tx,Iy) + kronecker(Ix,Ty)
    return(L)
}

circusTentQP <- function(nBase, nFolds ){
    tentbase <-makeTentEnsemble(nBase, nFolds )
    theta <- 1  # material coefficient on first order term in energy function; may want to adjust
    x <- tentbase$x
    y <- tentbase$y
    z <- tentbase$z
    Ny  <- nrow(z)
    Nx  <- ncol(z)
    N   <- Nx*Ny
    hx  <- 1/(Nx-1)
    hy  <- 1/(Ny-1)
    res <- list()
    res$Dmat <- hx*hy*Laplacian2D(1,1,Nx, Ny)                               
    res$dvec <- -theta*(hx*hy)*rep(1,N)       
    res$Amat <- t(diag(N))
    res$bvec <- matrix(z,N,1,byrow=FALSE)     # lower bound
    
    # QP data
    qp <- QP(Dmat= hx*hy*Laplacian2D(1,1,Nx, Ny), 
             dvec= -theta*(hx*hy)*rep(1,N), 
             Amat= t(diag(N)),
             bvec= matrix(z,N,1,byrow=FALSE)
    )
    return(qp)
}

