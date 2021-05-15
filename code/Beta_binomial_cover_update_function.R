beta_binomial_cover_update <- function(cover, presence_absence_data, n, remove_column) {
  library(fitdistrplus)
  # We remove the columns that are not wantend, this is done by only selecting
  # the columns from after the column number that the user want to remove. This is done for both
  # the cover data and the present/absence data
  cover_data <- cover[,(remove_column+1):ncol(cover)]
  freq_data <- presence_absence_data[,(remove_column+1):ncol(presence_absence_data)]

  # We now reomve the columns in cover data, that contain a specie, not present in any of the plots.
  # This is done by finding the species not present in any plot in the presence/absence data, and then
  # removing the corresponding columns from the cover data.
  for (species in colnames(cover_data)){
    if (sum(freq_data[,species]) == 0){
      cover = cover[,!(names(cover) %in% species)]
      cover_data <- cover_data[,!(names(cover_data) %in% species)]
    }

  }

  # We then create a data frame to hold the fitted parameter values for the prior distribution
  # for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))

  # We name the columns in the new data frame
  colnames(beta_fit) <- c("species","a", "b")

  # For each species we calculate the parameters for the fitted prior beta distribution and save them
  # in the data frame, beta_fit.
  for (name in colnames(cover_data)) {
    # We devide by n, to normalize the data, so each entri in the cover data, contains values between
    # 0 and 1.
    species <- cover_data[,name]/n

    # We remove all plots with 0 in frekvens for the given specie.
    beta_data <- species[freq_data[[name]] == 1]

    # We now want to fit a prior bete distribution for the given speice. We can only fit a beta
    # distribution if there is more than 1 unique value, and therefore we use the argument
    # length(unique(beta_data)) > 1. If there is not more than 1 unique value, we make the parameters
    # of the distribution 0. Else we use method og moment "mme" to fit a prior beta distribution.
    if (length(unique(beta_data)) > 1) {
      # Here we fit a beta distribution
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      # We then save the parameters in the data frame
      beta_fit[nrow(beta_fit) + 1,] <- c(name, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])

    }
    else {
      beta_fit[nrow(beta_fit) + 1,] <- c(name, 0,0)

    }
  }


  #For each cell in the input cover data we find the new abundance estimate and save it to this data frame.
  for (plot in 1:nrow(cover_data)) {

    #For each row (plot) we find the species that has a 1 in the presence/absence data.
    species_spotted_in_freq <- colnames(freq_data[c(freq_data[plot,]  == 1)])

    # We then find the species not in cover but is in the presence/absence data.
    not_in_cover <- setdiff(species_spotted_in_freq,colnames(cover_data))

    # We then choose the species that are both in cover and presence absence data.
    species_spotted_in_freq <- setdiff(species_spotted_in_freq, not_in_cover)


    #We calculate new abundance estimates for each spotted species in the plot
    for (species_spotted in species_spotted_in_freq ) {

      # We calculate the parameters of the posterior beta distribution for each specie in the plot
      alpha_post <- as.numeric((as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                                  as.numeric(cover_data[[species_spotted]][plot]) ))
      beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + n - as.numeric(cover_data[[species_spotted]][plot])

      # We find the mean value of the posteriro beta distribution for each specie in the plot,
      # with the parameters calculated above.
      mean_posterior <- (alpha_post)/(alpha_post+beta_post)

      # If the value is to small, we will change it to 0, else we time the value with n,
      # to unnormalise the value
      if(mean_posterior < 0.00001){
        cover[plot,species_spotted] <- 0
      }
      else{
        cover[plot,species_spotted] <- mean_posterior * n
      }

    }
  }
  #We then remove the columns from the cover data set, if the specie is not in any of the plots.
  for (species in colnames(cover_data)){
    if (sum(cover[,species]) == 0){
      cover = cover[,!(names(cover) %in% species)]
    }

  }

  # We return the updated cover data set which will be referenced to as beta binomial cover
  return(cover)
}
