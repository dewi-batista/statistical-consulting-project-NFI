library(readr)
library(dplyr)
library(mixtools)
library(evmix)

# Get Data and Names of markers and mixtures

data_mixtures = read.csv("C:\\Users\\rober\\Documents\\GitHub\\consulting_NFI\\mixtures.csv")
mixture_names = unique(data_mixtures$X)
markers_names = colnames(data_mixtures)[2:20]


gamma_mixtures = function(id_marker,id_mixture,max_components){
  print(mixture_names[id_mixture])
  print(markers_names[id_marker])
  mixture_marker_vector = data_mixtures[data_mixtures$X==mixture_names[id_mixture],][markers_names[id_marker]]
  # remove NA
  mixture_marker_vector = na.omit(mixture_marker_vector,as.numeric)
  # transform above dataframe to numeric 
  mixture_marker_vector = sapply(mixture_marker_vector,as.numeric)
  
  
  if (length(mixture_marker_vector)<5){
    return("Too little data")
  }
  
  
  # Perform model selection
  
  
  for(n_components in 1:max_components){
  
    output = capture.output({
    
      gamma_mixtures_results = gammamixEM(mixture_marker_vector,verb = FALSE,maxit=1e6,k=n_components,epsilon = 1e-08)
    })
  
    if (output[2] == "WARNING! NOT CONVERGENT! "){
      
    }

    
  
  }
}


#x = seq(1,10000,1)
#y = dmgamma(x , mgshape =gamma_results$gamma.pars[1,] , mgscale=gamma_results$gamma.pars[2,] , mgweight = gamma_results$lambda)
#hist(MUC4_semen_vaginal,prob = TRUE,breaks=100)
#lines(x,y)