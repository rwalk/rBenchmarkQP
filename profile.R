library(ggplot2)
library(reshape2)
source("./qp_examples.R")
source("./qp_solvers.R")

#######################################################
# Experiment setup and outputs
#######################################################
trial <- function(N, qp){
    # stage the solvers
    ipop_solve <- ipopStage(qp@Dmat, qp@dvec, qp@Amat, qp@bvec)
    quadprog_solve <- quadprogStage(qp@Dmat, qp@dvec, qp@Amat, qp@bvec)
    ipoptr_solve <- ipoptrStage(qp@Dmat, qp@dvec, qp@Amat, qp@bvec)
    osqp_solve <- osqpStage(qp@Dmat, qp@dvec, qp@Amat, qp@bvec)
    
    # timers
    quadprog.time <- system.time(quadprog.solution <- quadprog_solve())[3]
    ipop.time <- system.time(ipop.solution <- ipop_solve())[3]
    ipoptr.time <- system.time(ipoptr.solution <- ipoptr_solve())[3]
    osqp.time <- system.time(osqp.solution <- osqp_solve())[3]
    
    # error measurements
    ipop.diff <- sqrt(sum(quadprog.solution$solution - primal(ipop.solution))^2)
    ipoptr.diff <- sqrt(sum(quadprog.solution$solution - ipoptr.solution$solution)^2)
    osqp.diff <- sqrt(sum(quadprog.solution$solution - osqp.solution$x)^2)
    
    experiment <- c(N, 
                    quadprog.time, 
                    ipop.time, 
                    ipoptr.time, 
                    osqp.time, 
                    ipop.diff, 
                    ipoptr.diff, 
                    osqp.diff)
    
    names(experiment) <- c("problem.size", 
                           "quadprog.time", 
                           "ipop.time", 
                           "ipoptr.time", 
                           "osqp.time", 
                           "ipop.diff", 
                           "ipoptr.diff", 
                           "osqp.diff")
    return(experiment)
}

#######################################################
# Plotting
#######################################################
profilePlot <- function(times, filename, title){
    df <- melt(times, id = "problem.size", value.name = "time", variable.name="method")
    plt <- ggplot(df, aes(x=problem.size, y=time, color=method)) + 
        ggtitle(title)+
        geom_line(size=1.5) + 
        geom_point(size=5) + 
        scale_colour_manual(values=c("red","blue","yellow", "green")) +
        theme(legend.position="bottom", legend.title=element_blank())
    ggsave(filename, dpi=320)
}

# Plot a sparse matrix -- similar to MATLAB spy
spy <- function(A){
    M <- ifelse(A==0,0,1)
    image(t(M)[ncol(M):1,], col=c("white","blue"), axes=FALSE, pch = 3, lty=0)
    title("Sparse Matrix Plot")
    box()
}

######################################################
# Run the experiments
######################################################

#
# Circus tent example
#
circusTrial <- function(nFolds) trial((2**3*nFolds)**2, circusTentQP(2**3, nFolds))
circusExperiment <- as.data.frame(t(sapply(1:10, circusTrial)))

print(dput(circusExperiment))
print(circusExperiment)
write.csv(circusExperiment, file = "/tmp/qp_tent.data.csv", row.names = FALSE)
profilePlot(circusExperiment[,1:5], "/tmp/qp_tent.png", "Solve time vs problem size for circus tent problem")

#
# Random QP example
#
randomQPTrial <- function(N) trial(N, randomQP(N, 4))
problemSize <- seq(from = 500, to = 7500, by = 1000)
randomQPTrial <- as.data.frame(t(sapply(problemSize, randomQPTrial)))
print(dput(randomQPTrial))
print(randomQPTrial)
write.csv(randomQPTrial, file = "/tmp/random_qp.data.csv", row.names = FALSE)
profilePlot(randomQPTrial[,1:5], "/tmp/random_qp.png", "Solve time vs problem size for random dense QP")
