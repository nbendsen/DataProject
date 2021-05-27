
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



#Function for Species richness calculated using presence/absence data.
# We should only use columns that contains information on a specie, we can remove a number of columns,
# which does not contain informtaion on species, by setting remove_column to the number of columns,
# whished to be removed.
species_richness <- function(presence_absence_data, remove_column = NULL) {
  #Species reichness is calculated for each plot (row) and is just the sum of the row.

  if (is.null(remove_column)) {
    # We calulate the sum of each row. Each row represents a plot, so the sum will just be the number of species, present in the plot.
    species_richness <- rowSums(presence_absence_data)
  }
  else {
    # We calulate the sum of each row. Each row represents a plot, so the sum will just be the number of species, present in the plot.
    species_richness <- rowSums(presence_absence_data[, (remove_column+1):ncol(presence_absence_data)])
  }
  return(species_richness)
}


#Function for shannon index calculated using cover data.
# We should only use columns that contains information on a specie, we can remove a number of columns,
# which does not contain informtaion on species, by setting remove_column to the number of columns,
# whished to be removed.
shannon <- function(cover_data, remove_column = NULL){

  if(is.null(remove_column)){
    # If remove_column is NULL we use all og the columns
    data <- cover_data[,1:ncol(cover_data)]
  }
  else{
    # If remove_column is not NULL we remove the columns, not wanted
    data <- cover_data[,(remove_column+1):ncol(cover_data)]
  }

  # We calulate the sum of each row
  total <- rowSums(data)

  # We calculate the shannon index
  shannon_index <- rowSums(-data[,1:ncol(data)]/total *log(data[,1:ncol(data)]/total), na.rm = TRUE)
  return(shannon_index)
}

#Function for hill shonnon calculated using cover data.
# We should only use columns that contains information on a specie, we can remove a number of columns,
# which does not contain informtaion on species, by setting remove_column to the number of columns,
# whished to be removed.

hill_shannon <- function(cover_data, remove_column = NULL){
  if(is.null(remove_column)){

    # We use the shannon function to calculate the shannon indexes
    vector <- shannon(cover_data)
    # We calculate the hill shannon
    return(exp(vector))

  }
  else{
    # We use the shannon function to calculate the shannon indexes
    vector <- shannon(cover_data, remove_column = remove_column)

    # We calculate the hill shannon
    return(exp(vector))
  }


}

#Function for simpsons index calculated using cover data.
# We should only use columns that contains information on a specie, we can remove a number of columns,
# which does not contain informtaion on species, by setting remove_column to the number of columns,
# whished to be removed.
simpson <- function(cover_data, remove_column = NULL){

  if(is.null(remove_column)){
    # If remove_column is NULL we use all og the columns
    data <- cover_data[,1:ncol(cover_data)]
  }
  else{
    # If remove_column is not NULL we remove the columns, not wanted
    data <- cover_data[,(remove_column+1):ncol(cover_data)]
  }

  # We calulate the sum of each row
  total <- rowSums(data)

  # We calculate the simpson index
  simpson_index <- rowSums((data[,1:ncol(data)]/total)^2)

  simpson_index[is.nan(simpson_index)] <- 0

  return(simpson_index)
}


#Function for hill simpson calculated using cover data.
# We should only use columns that contains information on a specie, we can remove a number of columns,
# which does not contain information on species, by setting remove_column to the number of columns,
# whished to be removed.
hill_simpson <- function(cover_data, remove_column = NULL){
  if(is.null(remove_column)){
    # We use the simpson function to calculate the simpson indexes
    vector <- simpson(cover_data)

    # We then calculate the hill simpson
    v <- 1/vector
    v[is.infinite(v)] <- 0
    return(v)


  }
  else{
    # We use the simpson function to calculate the simpson indexes
    vector <- simpson(cover_data, remove_column = remove_column)

    # We then calculate the hill simpson
    a <- 1/vector
    a[is.infinite(a)] <- 0
    return(a)
  }


}


#A function to plot the different ways of measuring diversity against a value, i.e ph value,
# for each plot.
# description is the label for the x-axis
plot_diversity <- function(data, diversities, plot_info, description = NULL) {
  plot_data <-  gather(data, key = "type", value = "Diversity", diversities)
  ggplot(data = plot_data, mapping = aes(x = as.numeric(plot_data[[plot_info]]), y = Diversity)) +
    geom_point()+
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(vars(type),scales = "free_y") +

    #Label for the x-axis
    xlab(description) +

    # Label for the y-axis
    ylab("Diversity estimates") +
    theme_update(strip.placement = "outside", strip.text.x = element_text(face = "bold"), title = element_text(face = "bold"))
}



# Function for calculating the prior beta distributions for each specie.

beta <- function(cover, freq) {
  library(fitdistrplus)

  # We choose only the colums which contian information on species.
  cover_data <- cover[,4:ncol(cover)]
  freq_data <- freq[,4:ncol(freq)]

  # We remove the columns of the species in the cover data that is not
  # present in the presence/absence data
  for (ele in colnames(cover_data)){
    if (sum(freq_data[,ele]) == 0){
      cover = cover[,!(names(cover) %in% ele)]
      cover_data <- cover_data[,!(names(cover_data) %in% ele)]
    }

  }

  # Create data frame to hold the fitted values for the prior beta distribution for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))

  # Creating the column names for the new data frame
  colnames(beta_fit) <- c("species","a", "b")

  # For each species we calculate the parameters for the fitted prior beta distribution
  # and save them in the data frame.
  for (specie in colnames(cover_data)) {
    beta_data <- cover_data[,specie]/16

    # We remove all plots with 0 in frekvens for the given specie.
    beta_data <- beta_data[freq_data[[specie]] == 1]

    # We now want to fit a prior bete distribution for the given speice. We can only fit a beta
    # distribution if there is more than 1 unique value, and therefore we use the argument
    # length(unique(beta_data)) > 1. If there is not more than 1 unique value, we make the parameters
    # of the distribution 0. Else we use method og moment "mme" to fit a prior beta distribution.
    if (length(unique(beta_data)) > 1) {
      # Here we fit a beta distribution
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")

      # We then save the parameters in the data frame
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])

    }
    else {
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, 0,0)

    }


  }
  return(beta_fit)
}




different_diversities <- function(data_observed, data_new, plot, remove_column = NULL){
  observed <- data_observed[plot,(remove_column+1):ncol(data_observed)]
  update <- data_new[plot,(remove_column+1):ncol(data_new)]

  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("l", "Observed", "Updated")

  total_o <- rowSums(observed)
  total_u <- rowSums(update)

  for (l in seq(-1, -0.1, 0.1)){
    obs <- 0
    upt <- 0
    for (ele in 1:ncol(observed)){
      p <- observed[ele]/total_o
      if(p>0){
        r <- 1/p
        obs <- obs + p*(r^l)

      }
    }

    for (ele in 1:ncol(update)){
      p <- update[ele]/total_u
      if(p>0){
        r <- 1/p
        upt <- upt + p*(r^l)

      }
    }

    obs <- obs^(1/l)
    upt <- upt^(1/l)
    df[nrow(df)+1, ] <- c(l, obs, upt)

  }

  df[nrow(df)+1, ] <- c(0,hill_shannon(observed,0), hill_shannon(update,0))

  for (l in seq(0.1, 1, 0.1)){
    obs <- 0
    upt <- 0
    for (ele in 1:ncol(observed)){
      p <- observed[ele]/total_o
      if(p>0){
        r <- 1/p
        obs <- obs + p*(r^l)

      }
    }

    for (ele in 1:ncol(update)){
      p <- update[ele]/total_u
      if(p>0){
        r <- 1/p
        upt <- upt + p*(r^l)

      }
    }

    obs <- obs^(1/l)
    upt <- upt^(1/l)
    df[nrow(df)+1, ] <- c(l, obs, upt)

  }
  points <- df[df$l %in% c(-1,0,1),]
  ggplot(data = df, aes(x = l)) +
    geom_line(aes(y = Observed, colour = "Observed cover data"))+
    geom_line(aes(y = Updated, colour = "Beta binomial cover data"))+
    geom_point(data = points, mapping = aes(x = as.numeric(l), y = Updated),
               fill = "blue", shape=15,  size = 2, colour = "blue") +
    geom_point(data = points, mapping = aes(x = as.numeric(l), y = Observed),
               fill = "red", shape=17, size = 2, colour = "red") +
    scale_colour_manual("",
                        values = c("Observed cover data"="red",
                                   "Beta binomial cover data"="blue")) +

    labs(y = "Diversity", x = "Exponent l in Hill diversity formula")+
    ggtitle(sprintf("Comparison of diversity estimates for plot %d", plot))




}














