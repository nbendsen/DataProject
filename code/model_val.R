ppc <- function(m, freq_data, cover_data){
  library(fitdistrplus)
  
  
  cover_data <- cover[,4:ncol(cover)]
  freq_data <- freq[,4:ncol(freq)]
  
  
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))
  
  colnames(beta_fit) <- c("species","a", "b")
  
  for (specie in colnames(cover_data)) {
    
    beta_data <- cover_data[,specie]/16
    
    beta_data <- beta_data[freq_data[[specie]] == 1]
    
    
    if (length(unique(beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      alpha <- beta_data_fitted$estimate[1]
      beta <- beta_data_fitted$estimate[2]
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, alpha, beta)
    }
    else {
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, 0,0)
      
    }
  }
  
  
  n <- m
  
  species_spotted_in_frekvens <- colnames(freq_data[c(freq_data[n,]  == 1)])
  
  observed <- cover_data[n,c(species_spotted_in_frekvens)]
  tmp <- observed[observed > 0]
  T_static <- -sum(tmp/sum(observed) * log((tmp/sum(observed))))
  
  new_beta <- data.frame(matrix(ncol = 3, nrow = 0))
  
  colnames(new_beta) <- c("species","a", "b") 
  
  for (species_spotted in species_spotted_in_frekvens ) {
    
    alpha_post <- as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) + cover_data[[species_spotted]][n] 
    beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + 16 - cover_data[[species_spotted]][n]
    
    new_beta[nrow(new_beta) + 1,] <- c(species_spotted, alpha_post, beta_post)}
  
  
  shannon <- c()
  
  for (i in 1:1000){
    values <- c()
    for (ele in species_spotted_in_frekvens){
      a <- as.numeric(new_beta[new_beta$species == ele,]$a)
      b <- as.numeric(new_beta[new_beta$species ==ele,]$b)
      
      values <- c(values, rbeta(1,a,b))
      
      
    }
    tmp <- values[ values > 0.00001]
    total <- sum(tmp)
    shannon <- c(shannon,-sum(tmp/total * log((tmp/total)))) 
    
  }
  
  min_val <- min(T_static, min(shannon)) - 0.1
  max_val <- max(T_static, max(shannon)) + 0.1
  
  j <- length(shannon)
  pvalue <- 2*min(sum(shannon>= T_static)/j, sum(shannon<= T_static)/j)
  
  
  hist(shannon, xlim = c(min_val, max_val), main = sprintf("Histogram of simulated shannon indexes for plot %d", n), xlab = "Shannon indexes")
  legend("topright", legend = sprintf("Red line is \nobserved shannon index\n\n The posterior predictive \n p-value = %.3f", pvalue),box.lty=0)
  abline(v = T_static, col = "red" )

}