library(readr)
library(dplyr)
library(mixtools)
library(evmix)
library(stats)
set.seed(42) 

# Get Data and Names of markers and mixtures

# data_mixtures = read.csv("C:\\Users\\rober\\Documents\\GitHub\\consulting_NFI\\mixtures.csv")
data_mixtures = read.csv("./data/mixtures.csv")
mixture_names = unique(data_mixtures$X)
markers_names = colnames(data_mixtures)[2:16]


gamma_mixtures = function(id_mixture,id_marker,max_components){
  
  #
  # This function performs model selection in order to find the optimal number of components for a mixture of Gammas.
  # It will also plot an histogram of the true and generated overlapped with the fitted distribution. 
  #
  #
  # Args:
  #       id_mixture: integer. The index of the relevant fluid type in the list of fluids. The list:
  #                 ['Semen.fertile+Vaginal.mucosa' 'Saliva+Semen.fertile' 'Saliva+Vaginal.mucosa' 'Blood+Nasal.mucosa'
  #                  'Nasal.mucosa+Saliva''Vaginal.mucosa+Blood' 'Menstrual.secretion+Blood']
  #
  #       id_marker: integer. The index of the relevant marker in the list of markers from the data
  #                  The list: [HBB,ALAS2,CD93,HTN3,STATH,BPIFA1,MUC4,MYOZ1,CYP2B7P1,MMP10,MMP7,MMP11,SEMG1,KLK3,PRM1]
  #
  #       max_components: integer. The number of components to be used for the Gamma mixture model selection.
  # 
  # Return: best BIC
  
  
  
  mixture_marker_vector = data_mixtures[data_mixtures$X==mixture_names[id_mixture],][markers_names[id_marker]]
  # remove NA
  mixture_marker_vector = na.omit(mixture_marker_vector,as.numeric)
  # transform above dataframe to numeric 
  mixture_marker_vector = sapply(mixture_marker_vector,as.numeric)
  
  # if too little data is avaliable then running gamma mixtures does not produce any meaningful value
  # (particulalry when the whole vector is empty )
  if (length(mixture_marker_vector)<5){
    return(cat("Too little data for mixture ", mixture_names[id_mixture] , " and marker ", 
               markers_names[id_marker],"\n" ))
    
  }
  
  splitTrainTest = sample(c(TRUE, FALSE), length(mixture_marker_vector), replace=TRUE, prob=c(0.7,0.3))
  mixture_marker_vector_train = mixture_marker_vector[splitTrainTest]
  mixture_marker_vector_test = mixture_marker_vector[!splitTrainTest]
  
  #### MODEL SELECTION ###
  
  # initialize best parameters
  best_alpha = 0
  best_beta = 0
  best_lambda = 0
  best_BIC = 1e6 # very large BIC 
  
  for(n_components in 1:max_components){
    output = capture.output({
    
      gamma_mixtures_results = gammamixEM(mixture_marker_vector_train,verb = FALSE,maxit=1e6,k=n_components,epsilon = 1e-08)
      
    })
    
    
    #check if it converged
    if (isTRUE(output[2] == "WARNING! NOT CONVERGENT! ")==TRUE){
      print("no convergence")
      next
    }
    
    # check if BIC is better
    
    k = 3*n_components-1 # number of paramters to estimate (minus one because lambdas must sum to 1)
    
    current_BIC = k*log(length(mixture_marker_vector_train)) - 2*gamma_mixtures_results$loglik
    
    if(current_BIC<best_BIC){
      
      best_alpha = gamma_mixtures_results$gamma.pars[1,]
      best_beta = gamma_mixtures_results$gamma.pars[2,]
      best_lambda = gamma_mixtures_results$lambda
      best_BIC = current_BIC
      
    }
    
  }
  
  x = seq(min(mixture_marker_vector),max(mixture_marker_vector),1)
  y = dmgamma(x , mgshape =best_alpha , mgscale=best_beta , mgweight = best_lambda)
  new_data = rmgamma(length(mixture_marker_vector_test), mgshape = best_alpha, mgscale = best_beta,mgweight = best_lambda)
  
  par(mar = c(5, 4, 4, 8))
  
  hist(mixture_marker_vector,col=rgb(1,0,0,0.5),prob = TRUE,breaks=100,xlab = "Peak Values", 
       main=paste(mixture_names[id_mixture],",",markers_names[id_marker],"\nNumber of components = ", length(best_lambda)))
  
  hist(new_data,col=rgb(0,0,1,0.5),prob = TRUE,breaks=100,add=T)
  
  lines(x,y,col="green",lwd=4)
  
  legend("topright",xpd = TRUE, inset = c(-0.2, 0),bg = rgb(1, 1, 1, 0), legend=c("True Data","Generated Data","Gamma Mixture function"),
         col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),"green"), pt.cex=2, pch=15)
  
  return(best_BIC)
  
}



p_value_statistics = function(id_mixture,id_marker,n_components){
  
  
  #
  # This function generates 100 different datasets from a gamma mixture. Then performs the Kolomogrov-Smirnov test between the 
  # the test set and each new dataset. The corresponding p-values are saved in a vector. 
  #
  # Args:
  #       id_mixture: integer. The index of the relevant fluid type in the list of fluids. The list:
  #                 ['Semen.fertile+Vaginal.mucosa' 'Saliva+Semen.fertile' 'Saliva+Vaginal.mucosa' 'Blood+Nasal.mucosa'
  #                  'Nasal.mucosa+Saliva''Vaginal.mucosa+Blood' 'Menstrual.secretion+Blood']
  #
  #       id_marker: integer. The index of the relevant marker in the list of markers from the data
  #                  The list: [HBB,ALAS2,CD93,HTN3,STATH,BPIFA1,MUC4,MYOZ1,CYP2B7P1,MMP10,MMP7,MMP11,SEMG1,KLK3,PRM1]
  #
  #       n_components: integer. The number of components to be used for the Gamma mixture.
  # 
  # Return: median and minimum p-value
  
  
  p_values_vector = integer(100)

  mixture_marker_vector = data_mixtures[data_mixtures$X==mixture_names[id_mixture],][markers_names[id_marker]]
  # remove NA
  mixture_marker_vector = na.omit(mixture_marker_vector,as.numeric)
  # transform above dataframe to numeric 
  mixture_marker_vector = sapply(mixture_marker_vector,as.numeric)
  
  splitTrainTest = sample(c(TRUE, FALSE), length(mixture_marker_vector), replace=TRUE, prob=c(0.7,0.3))
  mixture_marker_vector_train = mixture_marker_vector[splitTrainTest]
  mixture_marker_vector_test = mixture_marker_vector[!splitTrainTest]
  gamma_mixtures_results = gammamixEM(mixture_marker_vector_train,verb = FALSE,maxit=1e6,k=n_components,epsilon = 1e-08)
  
  best_alpha = gamma_mixtures_results$gamma.pars[1,]
  best_beta = gamma_mixtures_results$gamma.pars[2,]
  best_lambda = gamma_mixtures_results$lambda
  
  for(i in 1:100){
    
    new_data = rmgamma(length(mixture_marker_vector_test), mgshape = best_alpha, mgscale = best_beta,mgweight = best_lambda)
    ks_test_new_data = ks.test(new_data,mixture_marker_vector_test)
    p_values_vector[i] = ks_test_new_data$p.value
  }
  
  
  return(c(median(p_values_vector),min(p_values_vector)))
  
}











