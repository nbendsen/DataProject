ppc <- function(m, freq_data, cover_data){
  library(fitdistrplus)

  #We read the cover data and the presence/absence data without the first 4 columns,
  # as they do not contains information on species
  cover_data <- cover[,4:ncol(cover)]
  freq_data <- freq[,4:ncol(freq)]
  #We make a dataframa for the parameters of the prior distribution for each plot,
  #so it is possible to save the parameters, and not calculate them when making the posterior for each plot.
  #Each row will contain the number/name of the specie and its corresponding parameters for the prior

  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))
  # We name the columns
  colnames(beta_fit) <- c("species","a", "b")
  #Here we calculate the parameters for the prior distributions of each specie:
  for (specie in colnames(cover_data)) {
    #First we normalize. Since there is in total used 16 pins for each plot,
    #we will divide the entries in the cover data by 16
    beta_data <- cover_data[,specie]/16

    #Now we remove the plots where the specie is not present.
    #This can be done by using the information from the presence, absence data.
    #If it contains a 1, then the specie is present in the plot, if 0 it is absent.
    beta_data <- beta_data[freq_data[[specie]] == 1]

    #If the specie is not present in any of the plots,
    #we do not have information to make a prior distribution for it,
    #and will just give it parameters a=0 and b=0 as seen in the else clause.

    if (length(unique(beta_data)) > 1) {
      #We use the method of moments to fit the prior beta distribution
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      #The parameters are added to the dataframe
      alpha <- beta_data_fitted$estimate[1]
      beta <- beta_data_fitted$estimate[2]
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, alpha, beta)
    }
    else {
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, 0,0)

    }
  }

  #n is the row number of the plot we are working with
  n <- m
  #We define which species are in present in the plot
  species_spotted_in_frekvens <- colnames(freq_data[c(freq_data[n,]  == 1)])
  # we remove the columns, that are not representing species
  observed <- cover_data[n,c(species_spotted_in_frekvens)]
  tmp <- observed[observed > 0]
  T_static <- exp(-sum(tmp/sum(observed) * log((tmp/sum(observed)))))
  #We make a dataframe to save the parameters of the posterior for each spotted specie in the plot
  new_beta <- data.frame(matrix(ncol = 3, nrow = 0))

  colnames(new_beta) <- c("species","a", "b")

  for (species_spotted in species_spotted_in_frekvens ) {
    #We define the parameters for the posterior
    alpha_post <- as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) + cover_data[[species_spotted]][n]
    beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + 16 - cover_data[[species_spotted]][n]

    #We add the parameters to the dataframe
    new_beta[nrow(new_beta) + 1,] <- c(species_spotted, alpha_post, beta_post)}

  #We make a vector, to save the shannon indexes produces in each iteration
  shannon <- c()

  for (i in 1:1000){
    #Vector for saving the random generated values from the posterior of each specie
    values <- c()
    for (ele in species_spotted_in_frekvens){
      #These are the parameters, for the posterior of the specie
      a <- as.numeric(new_beta[new_beta$species == ele,]$a)
      b <- as.numeric(new_beta[new_beta$species ==ele,]$b)
      # We draw a random number from a beta distribution with the parameters
      # for that specie and add it to the vector
      values <- c(values, rbeta(1,a,b))


    }
    #We remove the values that are to small
    tmp <- values[ values > 0.00001]
    total <- sum(tmp)
    #We calculate the Hill Shannon diversity
    shannon <- c(shannon, exp(-sum(tmp/total * log((tmp/total)))))

  }

  min_val <- min(T_static, min(shannon)) - 0.1
  max_val <- max(T_static, max(shannon)) + 0.1

  j <- length(shannon)
  pvalue <- 2*min(sum(shannon>= T_static)/j, sum(shannon<= T_static)/j)

  #This is the code that produces the histogram over the generated data
  hist(shannon, xlim = c(min_val, max_val), breaks = 15, main = sprintf("Histogram of simulated Hill Shannon diversity for plot %d", n), xlab = "Hill Shannon diversity")
  legend("topright", legend = sprintf("Red line is observed \nHill Shannon diversity \n\n The posterior predictive \n p-value = %.3f", pvalue),box.lty=0)
  abline(v = T_static, col = "red" )

}
