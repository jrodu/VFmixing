#' Generate Individual Cubic Trajectories
#'
#' @param startTime a numeric scalar. The center time for the generating cubic function.
#' @param endTime a numeric scalar. The ending time for the generating cubic function. 
#' If another cubic were to be generated after, this would be the starting time for the next.
#' @param prevTraj an nx4 numeric matrix. Each row corresponds to the coefficients of the previous cubic functions. 
#' This argument is not needed if initTraj = TRUE. 
#' @param prevStartTime a numeric scalar. The centering time for the previous cubic function.
#' This argument is not needed if initTraj = TRUE.
#' @param initTraj a boolean. TRUE if the generating cubic is the first, FALSE if there are cubic functions preceding. 
#' @param beta4 a numeric vector. The constant coefficient for the initial cubic function. 
#' This argument is not needed if initTraj = TRUE.
#' @param endSD a non-negative numeric scalar. 
#' The data generating process says the distribution on the difference between the starting and end value of the generated trajectory is normal with mean 0 and standard deviation endSD. 
#'
#' @returns an nx4 matrix. The generated cubic functions with each row corresponding to the coefficients. 
#' @export
#'
#' @examples
#' 
#' startTime <- 2
#' endTime <- 3
#' prevTraj <- matrix(c(0,0,0,1), nrow = 1, byrow = TRUE)
#' prevStartTime <- 1
#' initTraj <- FALSE
#' endSD <- 0.5
#' 
#' generateIndividualTrajectory(startTime = startTime, endTime = endTime, 
#' prevTraj = prevTraj, prevStartTime = prevStartTime, initTraj = initTraj, endSD = endSD)
#' 
#' 
generateIndividualTrajectory = function(startTime, endTime, prevTraj, prevStartTime, initTraj = F, beta4, endSD = 0.5){
  
  curTimeDiff = endTime - startTime
  
  if(initTraj){
    mu = c(0,0,0)
    sigma = (2*endSD^2)*matrix(c((curTimeDiff^2+2)/(curTimeDiff^6)        , (-3*curTimeDiff^2 - 6)/(2*curTimeDiff^5), (1)/(2*curTimeDiff^2),
                                 (-3*curTimeDiff^2 - 6)/(2*curTimeDiff^5) , (5*curTimeDiff^2 + 9)/(2*curTimeDiff^4) , (-1)/(curTimeDiff),
                                 (1)/(2*curTimeDiff^2)                    , (-1)/(curTimeDiff)                      , 1/2), 
                               nrow = 3, byrow = T)
    
    sampTraj = MASS::mvrnorm(n = length(beta4), mu = mu, Sigma = sigma)
    curTraj = cbind(matrix(sampTraj, nrow = length(beta4)), beta4,deparse.level = 0)
    

  } else{
    
    prevTimeDiff = startTime - prevStartTime
    
    squig = t(t(3*prevTraj[,1]*prevTimeDiff^2 + 2*prevTraj[,2]*prevTimeDiff + prevTraj[,3]))
    
    mu = cbind((squig)/(curTimeDiff^2), (-2*squig)/(curTimeDiff))
    sigma = (2*endSD^2)*matrix(c((curTimeDiff^2+4)/(2*curTimeDiff^6) , (-curTimeDiff^2-6)/(2*curTimeDiff^5),
                                 (-curTimeDiff^2-6)/(2*curTimeDiff^5), (curTimeDiff^2+9)/(2*curTimeDiff^4)), 
                               nrow = 2, byrow = T)
    
    sampTraj = MASS::mvrnorm(n = nrow(prevTraj),mu = rep(0,2), Sigma = sigma) + mu
    
    curTraj = generateC1Smoothness(prevCubic = prevTraj, curCubic = matrix(sampTraj, ncol = 2), 
                                  prevStartTime = prevStartTime, curStartTime = startTime, 
                                  meetTime = startTime)
    
  }
  
  colnames(curTraj) = c("Beta1","Beta2","Beta3","Beta4")
  curTraj
  
}




#' Generate Full Trajectory Function
#'
#' @param TimeSplits a numeric vector. The spline node locations with the starting and ending times. 
#' @param numModels a numeric scalar. The number of trajectories to generate.
#' @param sigma a non-negative numeric scalar. The standard deviation of the normal distribution generating the starting values of the trajectories. Not needed if startWeights is specified.
#' @param k a non-negative numeric scalar. The standard deviation of the induced normal distribution on difference between the starting and ending values of the individual cubic functions.
#' This is the equivalent of endSD in generateIndividualTrajectory function.
#' @param startWeights a numeric vector. The starting values of the trajectory functions. 
#'
#' @returns a list. Each element is a coefficient matrix corresponding to a segment of the spline. 
#' @export
#'
#' @examples
#' 
#' generateFullTrajectory(TimeSplits = 0:5, numModels = 2, sigma = 0.1, k = 0.5)
#' 
generateFullTrajectory = function(TimeSplits, numModels, sigma, k, startWeights){
  if(missing(startWeights)){
    startWeights = stats::rnorm(numModels, mean = 1/numModels, sd = sigma)
  }
  
  randTrajList = list()
  
  randTrajList[[1]] = generateIndividualTrajectory(startTime = TimeSplits[1], endTime = TimeSplits[2], 
                                                   endSD = k, initTraj = T, prevTraj = NULL, prevStartTime = NULL,beta4 = startWeights)
  
  if(length(TimeSplits)>2){
    
    for(i in 2:(length(TimeSplits)-1)){
      randTrajList[[i]] = generateIndividualTrajectory(startTime = TimeSplits[i], endTime = TimeSplits[i+1], endSD = k, initTraj = F, prevStartTime = TimeSplits[i-1],prevTraj = randTrajList[[i-1]])
    }
    
  }
  
  randTrajList
  
}


