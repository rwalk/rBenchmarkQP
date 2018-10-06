require(quadprog)
require(kernlab)
require(ipoptr)
require(osqp)

######################################################
# QP solver staging
######################################################

# quadprog
quadprogStage <- function(Dmat, dvec, Amat, bvec){
    f <- function(){
        # add posiitivy constraint
        N <- length(dvec)
        Amat <- rbind(Amat, diag(1,N,N))
        bvec <- rbind(bvec, matrix(0,N,1))
        
        # note, quadprog constraint set requires a transpose
        sol <- solve.QP(Dmat, dvec, t(Amat), bvec)
        return(sol)
    }
    return(f)
}

# ipop
ipopStage <- function(Dmat, dvec, Amat, bvec){ 
    # solve the QP with kernlab solver
    # ipop solves the quadratic programming problem :
    # \min(c'*x + 1/2 * x' * H * x)
    # subject to: 
    # b <= A * x <= b + r
    # l <= x <= u
    ub <- 1e8
    n <- length(dvec)
    cc <- -dvec     # note the sign difference relative to quadprog
    H <- Dmat
    b <- bvec
    r <- matrix(ub, nrow=length(bvec), ncol=1)
    A <- Amat  
    
    # box constraint over x
    l <- matrix(0, nrow=n, ncol=1)
    u <- matrix(ub, nrow=n, ncol=1)
    
    f <- function(){
        # note, solver seems to be pretty sensitive to the margin size
        # if you see errors about conditioning try tweaking the margin size
        # up or down.
        sol <- ipop(cc, H, A, b, l, u, r, margin=1e-4)
        return(sol)
    }
    return(f)
}

# ipoptr
ipoptrStage <- function(Dmat, dvec, Amat, bvec, ub=1e10){
    # Stage the quadratic program
    #
    # min -d^T x + 1/2 x^T D x
    #   s.t. t(A)%*%x>= b
    #
    # for ipoptr given the matrix inputs.
    
    n <- length(bvec)
    N <- length(dvec)
    
    # Jacobian structural components
    j_vals <- unlist(Amat[Amat!=0])
    j_mask <- make.sparse(ifelse(Amat!=0, 1, 0))
    
    # Hessian structural components
    h_mask <- make.sparse(ifelse(upper.tri(Dmat, diag=TRUE) & Dmat!=0, 1, 0))
    h_vals <- do.call(c, sapply(1:length(h_mask), function(i) Dmat[i,h_mask[[i]]]))
    
    # build the ipoptr inputs
    eval_f <- function(x) return(-t(dvec)%*%x + 0.5*t(x)%*%Dmat%*%x)
    eval_grad_f <- function(x) return(-dvec + Dmat%*%x)
    eval_h <- function(x, obj_factor, hessian_lambda) return(obj_factor*h_vals)
    eval_h_structure <- h_mask
    eval_g <- function(x) return(Amat%*%x)
    eval_jac_g <- function(x) return(j_vals)
    eval_jac_g_structure <- j_mask
    constraint_lb <- bvec
    constraint_ub <- rep( ub, n)
    lb <- rep(0, N)
    
    # wrapper function that will solves the ipoptr problem when called
    f <- function() {
        # initialize with the global unconstrained minimum; this is a reasonable guess and consistent with quadprog.
        # NOTE: This will only work if lb <= x0 <= ub.  If this is not the case, 
        # use x0 = lb can be used instead.
        x0 <- solve(Dmat, dvec)
        if(!all(x0>=lb)) lb=x0
            
        # call the solver
        res <- ipoptr(x0 = x0, 
                      eval_f = eval_f,
                      eval_grad_f = eval_grad_f,
                      lb = lb,
                      eval_g = eval_g, 
                      eval_jac_g = eval_jac_g,
                      eval_jac_g_structure = eval_jac_g_structure,
                      constraint_lb = constraint_lb,
                      constraint_ub = constraint_ub,
                      eval_h = eval_h,
                      eval_h_structure = eval_h_structure)
        return(res)
    }
    return(f)
}

# OSQP
osqpStage <- function(Dmat, dvec, Amat, bvec){
    f <- function(){
        # add posiitivy constraint
        N <- length(dvec)
        Amat <- rbind(Amat, diag(1,N,N))
        bvec <- rbind(bvec, matrix(0,N,1))
        settings <- osqpSettings(verbose = FALSE, eps_abs=1e-8, eps_rel = 1e-8)
        sol <- solve_osqp(Dmat, -dvec, Amat, bvec, pars=settings)
        return(sol)
    }
    return(f)
}

