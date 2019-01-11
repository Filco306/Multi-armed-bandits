####################################
######## Initialize functions ######
####################################

Get_Convergence = function(mus) {
  sortedMaxRewards = sort(mus, decreasing = TRUE)
  
  return(sortedMaxRewards[1] - sortedMaxRewards[2])
}
normalize = function(vec) {
  
  print(sd(vec))
  if (sd(vec) != 0) {
    return((vec - mean(vec)/sd(vec)))
  } else {
    return(vec)
  }
  
}
Convergence_Means <- function (samples) {
  means <- numeric()
  for (i in 1:length(samples)) {
    means[i] <- mean(samples[1:i])
  }
  return(means)
}

ret_epsilon = function(mode = 1, t = NULL, c = 1) {
  if (mode == 1) {
    return(0.5)
  } else if (mode == 2) {
    return(min(1,c/sqrt(t)))
  } else if (mode == 3) {
    return(min(1,c/log(t)))
  } else if (mode == 4) {
    # Another suggestion?
  }
}

epsilonGreedy = function(k, mus, sigma2) {
  
  epsilon = ret_epsilon(mode = 1)
  
  # Set all mu estimates to 0.5 first.
  muEstimates = rep(0.5, k)
  
  # Set probability to pick each arm to epsilon/k, 
  probabilities = rep(epsilon/k, k)
  
  # Set the best arm's probability to 1 - epsilon + e/k, i.e. add 1 - epsilon
  # In this case, all
  # Initialize vectors to store information
  rewards = as.numeric()
  armChosen = as.numeric()
  t = 1
  while (t <= T) {
    
    # update epsilon 
    epsilon = ret_epsilon(epsilonMode, t = t)
    
    
    # Sample an arm
    armChosen[t] = sample(seq(1, k, 1), size = 1, prob = probabilities)
    
    rewards[t] = rnorm(n=1,mean = mus[armChosen[t]], sd = sqrt(sigma2))
    
    
    # Update muEstimates
    muEstimates[armChosen[t]] = sum(rewards[armChosen == armChosen[t]])/sum(armChosen == armChosen[t])
    
    # Update probabilities
    probabilities = rep(epsilon/k, k)
    # Set the best arm's probability to 1 - epsilon + e/k, i.e. add 1 - epsilon
    probabilities[which.max(muEstimates)] = probabilities[which.max(muEstimates)] + 1 - epsilon
    
    t = t + 1
  }
  
  data = data.frame(armChosen, rewards)
  colnames(data) = c("ArmChosen", "Rewards")
  return(data)
}

SoftMax = function(k, mus, sigma2, c = NULL) {
  muEstimates = rep(0.5, k)
  
  # Initialize vectors to store information
  rewards = as.numeric()
  armChosen = as.numeric()
  n_k_t = rep(0,k)
  
  for (t in 1:T) {
    if (!is.null(c)) {
      eta_t = c/t
    } else {
      eta_t = 0.9
    }
    
    
    probs = exp(eta_t*muEstimates)
    
    
    if(any(is.na(probs))) {
      print(probs)
      print(muEstimates)
      muEstimates = normalize(muEstimates)
      print(muEstimates)
      probs = exp(eta_t*muEstimates)/sum(exp(eta_t*muEstimates))
      print(probs)
    }
    armChosen[t] = sample(seq(1,k,1), size = 1, prob = probs)
    rewards[t] = rnorm(n=1,mean = mus[armChosen[t]], sd = sqrt(sigma2))
    n_k_t[armChosen[t]] = n_k_t[armChosen[t]] + 1
    muEstimates[armChosen[t]] = sum(rewards[armChosen == armChosen[t]])/(n_k_t[armChosen[t]] + max(0, 1 - n_k_t[armChosen[t]]))
  }
  data = data.frame(armChosen, rewards)
  colnames(data) = c("ArmChosen", "Rewards")
  return(data)
}


UCB = function(k, mus, sigma2) {
  
  muEstimates = rep(0.5, k)
  
  # Initialize vectors to store information
  rewards = as.numeric()
  armChosen = as.numeric()
  n_k_t = rep(0,k)
  
  t = 1
  while (t <= T) {
    
    # Generate the ucbs. We need the indices for this, hence the for loop.
    ucbs = as.numeric()
    for (i in 1:k) {
      
      ucbs[i] = muEstimates[i] + sqrt(2*log(t)/(n_k_t[i] + max(0, 1 - n_k_t[i])))
    }

    # Choose arm from the one with highest UCB
    armChosen[t] = which.max(ucbs)
    # update n_k_t
    n_k_t[armChosen[t]] = n_k_t[armChosen[t]] + 1
    
    # Sample reward
    rewards[t] = rnorm(n=1, mean = mus[armChosen[t]], sd = sqrt(sigma2))
    
    
    # Update muEstimates
    
    muEstimates[armChosen[t]] = sum(rewards[armChosen == armChosen[t]])/(n_k_t[armChosen[t]] + max(0, 1 - n_k_t[armChosen[t]]))
    
    t = t + 1
    
  }
  
  data = data.frame(armChosen, rewards)
  colnames(data) = c("ArmChosen", "Rewards")
  return(data)
}

calc_regret = function(muStar, armChosen, mus) {
  
  return(apply(as.matrix(armChosen), 2, function(row, muStar) {
    return(muStar - mus[armChosen])
  }, muStar))
}

cumRegret = function(regrets) {
  cumRegrets = as.numeric()
  for (i in 1:length(regrets)) {
    cumRegrets[i] = sum(regrets[1:i])
  }
  return(cumRegrets)
}

analyze = function(k, sigma2, regretsEps1, regretsEps2, regretsEps3, regretsUCB, regretsSM1, regretsSM10, avgBestArmEps1, avgBestArmEps2, avgBestArmEps3, avgBestArmUCB, avgBestArmSM1, avgBestArmSM10) {
  xSeq = seq(1, length(regretsEps1),1)
  plot(x=xSeq, y = Convergence_Means(regretsEps1), type = "l", lwd = 2, xlab = "Iteraciones", ylab = "Arrepentimiento promedio", col = "blue", ylim = c(0,max(Convergence_Means(regretsEps1), Convergence_Means(regretsUCB))),main = paste("K = ",k, ", ", expression(sigma), " = ", sigma2))
  lines(Convergence_Means(regretsUCB), lwd = 2, col = "red")
  lines(Convergence_Means(regretsEps2), lwd = 2, col = "green")
  lines(Convergence_Means(regretsEps3), lwd = 2, col = "black")
  lines(Convergence_Means(regretsSM1), lwd = 2, col = "purple")
  lines(Convergence_Means(regretsSM10), lwd = 2, col = "yellow")
  cumulativeRegretEpsilon1 = cumRegret(regretsEps1)
  cumulativeRegretEpsilon2 = cumRegret(regretsEps2)
  cumulativeRegretEpsilon3 = cumRegret(regretsEps3)
  cumulativeRegretUCV = cumRegret(regretsUCB)
  cumulativeRegretSM1 = cumRegret(regretsSM1)
  cumulativeRegretSM10 = cumRegret(regretsSM10)
  
  legend("top",
         col = c("blue", "red", "green", "black", "purple", "yellow"), 
         lwd = c(2,2), 
         legend = c(expression(paste(epsilon,"-greedy, ", epsilon, " = 0.5")),
                    "UCB",
                    expression(paste(epsilon," = 1/",sqrt(t))),
                    expression(paste(epsilon, "= 1/log(t)")),
                    expression(paste("Softmax, ",eta," = 0.9")),
                    expression(paste("Softmax, ",eta," = 10/t"))
         ))
  plot(x=xSeq, y = cumulativeRegretEpsilon1, type = "l", lwd = 2, xlab = "Iteraciones", ylab = "Arrepentimiento cumulativo", col = "blue", ylim = c(min(cumulativeRegretUCV, cumulativeRegretEpsilon1),max(cumulativeRegretUCV, cumulativeRegretEpsilon1)),main = paste("K = ",k, ", ", expression(sigma), " = ", sigma2))
  lines(cumulativeRegretUCV, lwd = 2, col = "red")
  lines(cumulativeRegretEpsilon2, lwd = 2, col = "green")
  lines(cumulativeRegretEpsilon3, lwd = 2, col = "black")
  lines(cumulativeRegretSM1, lwd = 2, col = "purple")
  lines(cumulativeRegretSM10, lwd = 2, col = "yellow")
  legend("top",
         col = c("blue", "red", "green", "black", "purple", "yellow"), 
         lwd = c(2,2), 
         legend = c(expression(paste(epsilon,"-greedy, ", epsilon, " = 0.5")),
                    "UCB",
                    expression(paste(epsilon," = 1/",sqrt(t))),
                    expression(paste(epsilon, "= 1/log(t)")),
                    expression(paste("Softmax, ",eta," = 0.9")),
                    expression(paste("Softmax, ",eta," = 10/t"))
                    ))
  
  ## Plot best arm
  
  plot(x=xSeq, y = avgBestArmEps1, type = "l", lwd = 2, xlab = "Iteraciones", ylab = "Veces mejor brazo elegido", col = "blue", ylim = c(0,max(avgBestArmEps1, avgBestArmEps2, avgBestArmEps3, avgBestArmUCB)),main = paste("K = ",k, ", ", expression(sigma), " = ", sigma2))
  lines(avgBestArmUCB, lwd = 2, col = "red")
  lines(avgBestArmEps2, lwd = 2, col = "green")
  lines(avgBestArmEps3, lwd = 2, col = "black")
  lines(avgBestArmSM1, lwd = 2, col = "purple")
  lines(avgBestArmSM10, lwd = 2, col = "yellow")
  legend("top",
         col = c("blue", "red", "green", "black", "purple", "yellow"), 
         lwd = c(2,2), 
         legend = c(expression(paste(epsilon,"-greedy, ", epsilon, " = 0.5")),
                    "UCB",
                    expression(paste(epsilon," = 1/",sqrt(t))),
                    expression(paste(epsilon, "= 1/log(t)")),
                    expression(paste("Softmax, ",eta," = 0.9")),
                    expression(paste("Softmax, ",eta," = 10/t"))
         ))
}

bestArmPlayed = function(armPlayed, mus) {
  percArmPlayed = as.numeric()
  for (i in 1:length(armPlayed)) {
    tmp = armPlayed[1:i]
    percArmPlayed[i] = sum(tmp == which.max(mus))
  }
  return(percArmPlayed)
}

############################################
######## Initialize and run ################
############################################

set.seed(123)
# Set parameters

K = c(5,10,25,50)

sigma2 = c(0.01^2,0.1^2,1^2, 10^2)

## Set algorithm specifications
# Set epsilon mode. Sets which way the epsilon is decided. 1 for fixed, 2 for epsilon = 1/t, 3 for epsilon = 1/sqrt(t), 4 for epsilon = 1/log(t)
epsilonMode = 1
c = 1
#Setting number of max iterations and convergence criteria. Problem will run until convergence or for the max nr iterations specified.
T = 2000
nRuns = 1
convergence = Get_Convergence(mus)
# Generate results
for (k in K) {
  for (sig in sigma2) {
    for (i in 1:nRuns) {
      # Get results for the run in terms of ArmChosen and reward
      mus = runif(k, min = 0, max = 1)
      epsilonMode = 1
      epsilonResults1 = epsilonGreedy(k, mus, sig)
      epsilonMode = 2
      epsilonResults2 = epsilonGreedy(k, mus, sig)
      epsilonMode = 3
      epsilonResults3 = epsilonGreedy(k, mus, sig)
      UCBResults = UCB(k, mus, sig)
      if (sig != 10^2) {
        softMaxC1 = SoftMax(k, mus, sig)
        softMaxC10 = SoftMax(k, mus, sig, c = 10)
      }
      # Now, convert to AVERAGE regret over the 1000 iterations. 
      if (i == 1) {
        avgRegretEps1 = calc_regret(mus[which.max(mus)], epsilonResults1$ArmChosen, mus)
        avgRegretEps2 = calc_regret(mus[which.max(mus)], epsilonResults2$ArmChosen, mus)
        avgRegretEps3 = calc_regret(mus[which.max(mus)], epsilonResults3$ArmChosen, mus)
        avgRegretUCB = calc_regret(mus[which.max(mus)], UCBResults$ArmChosen, mus)
        avgRegretSM1 = calc_regret(mus[which.max(mus)], softMaxC1$ArmChosen, mus)
        avgRegretSM10 = calc_regret(mus[which.max(mus)], softMaxC10$ArmChosen, mus)
        
        avgBestArmEps1 = bestArmPlayed(epsilonResults1$ArmChosen, mus)
        avgBestArmEps2 = bestArmPlayed(epsilonResults2$ArmChosen, mus)
        avgBestArmEps3 = bestArmPlayed(epsilonResults3$ArmChosen, mus)
        avgBestArmUCB = bestArmPlayed(UCBResults$ArmChosen, mus)
        avgBestArmSM1 = bestArmPlayed(softMaxC1$ArmChosen, mus)
        avgBestArmSM10 = bestArmPlayed(softMaxC10$ArmChosen, mus)
      } else {
        avgRegretEps1 = avgRegretEps1*(i-1)/i + calc_regret(mus[which.max(mus)], epsilonResults1$ArmChosen, mus)/i
        avgRegretEps2 = avgRegretEps2*(i-1)/i + calc_regret(mus[which.max(mus)], epsilonResults2$ArmChosen, mus)/i
        avgRegretEps3 = avgRegretEps3*(i-1)/i + calc_regret(mus[which.max(mus)], epsilonResults3$ArmChosen, mus)/i
        avgRegretUCB = avgRegretUCB*(i-1)/i + calc_regret(mus[which.max(mus)], UCBResults$ArmChosen, mus)/i
        avgRegretSM1 = avgRegretSM1*(i-1)/i + calc_regret(mus[which.max(mus)], softMaxC1$ArmChosen, mus)/i
        avgRegretSM10 = avgRegretSM10*(i-1)/i + calc_regret(mus[which.max(mus)], softMaxC10$ArmChosen, mus)/i
        
        avgBestArmEps1 = avgBestArmEps1*(i-1)/i + bestArmPlayed(epsilonResults1$ArmChosen,mus)/i
        avgBestArmEps2 = avgBestArmEps2*(i-1)/i + bestArmPlayed(epsilonResults2$ArmChosen,mus)/i
        avgBestArmEps3 = avgBestArmEps3*(i-1)/i + bestArmPlayed(epsilonResults3$ArmChosen,mus)/i
        avgBestArmUCB = avgBestArmUCB*(i-1)/i + bestArmPlayed(UCBResults$ArmChosen,mus)/i
        avgBestArmSM1 = avgBestArmSM1*(i-1)/i + bestArmPlayed(softMaxC1$ArmChosen,mus)/i
        avgBestArmSM10 = avgBestArmSM10*(i-1)/i + bestArmPlayed(softMaxC10$ArmChosen,mus)/i
      }
      
      print(i)
    }
    
    analyze(k, sqrt(sig), avgRegretEps1, avgRegretEps2, avgRegretEps3, avgRegretUCB, avgRegretSM1, avgRegretSM10, avgBestArmEps1, avgBestArmEps2, avgBestArmEps3, avgBestArmUCB, avgBestArmSM1, avgBestArmSM10)
    
  }
}
