bayesian_method <- function(cover, presence_absence_data, n, remove_column) {
  library(fitdistrplus)
  cover_data <- cover[,(remove_column+1):ncol(cover)]
  freq_data <- presence_absence_data[,(remove_column+1):ncol(presence_absence_data)]


  for (species in colnames(cover_data)){
    if (sum(freq_data[,species]) == 0){
      cover = cover[,!(names(cover) %in% species)]
      cover_data <- cover_data[,!(names(cover_data) %in% species)]
    }

  }

  #create data frame to hold the fitted values for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))

  colnames(beta_fit) <- c("species","a", "b")

  # for Each species calculate the shape parameter for the fitted beta distribution and save them in a data frame.
  for (name in colnames(cover_data)) {
    species <- cover_data[,name]/n

    #remove all plots with 0 in frekvens.
    beta_data <- species[freq_data[[name]] == 1]


    if (length(unique(beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      beta_fit[nrow(beta_fit) + 1,] <- c(name, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])

    }
    else {
      beta_fit[nrow(beta_fit) + 1,] <- c(name, 0,0)

    }
  }


  #For each cell in the input cover data we find the new abundance estimate and save it to this data frame
  for (plot in 1:nrow(cover_data)) {

    #For each row we find the species that has a 1 in the presence/absence data
    species_spotted_in_freq <- colnames(freq_data[c(freq_data[plot,]  == 1)])

    not_in_cover <- setdiff(species_spotted_in_freq,colnames(cover_data))

    species_spotted_in_freq <- setdiff(species_spotted_in_freq, not_in_cover)


    #We calculate new abundance estimate for each spotted species in the plot
    for (species_spotted in species_spotted_in_freq ) {


      alpha_post <- as.numeric((as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                                  as.numeric(cover_data[[species_spotted]][plot]) ))
      beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + n - as.numeric(cover_data[[species_spotted]][plot])

      mean_posterior <- (alpha_post)/(alpha_post+beta_post)

      if(mean_posterior < 0.00001){
        cover[plot,species_spotted] <- 0
      }
      else{
        cover[plot,species_spotted] <- mean_posterior * n
      }

    }
  }

  for (species in colnames(cover_data)){
    if (sum(cover[,species]) == 0){
      cover = cover[,!(names(cover) %in% species)]
    }

  }


  return(cover)
}



#Function for Species richness from frequency data.
#the species indicator is a substring that is present in only the the columns
# with data for species. i.e not plot, site or year. If species indicator is left empty, all columns will be treated as species columns.

species_richness <- function(presence_absence_data, remove_column = NULL) {

  if (is.null(remove_column)) {
    species_richness <- rowSums(presence_absence_data)
  }
  else {
    species_richness <- rowSums(presence_absence_data[, (remove_column+1):ncol(presence_absence_data)])
  }
  return(species_richness)
}



shannon <- function(cover_data, remove_column = NULL){

  if(is.null(remove_column)){
    data <- cover_data[,1:ncol(cover_data)]
  }
  else{
    data <- cover_data[,(remove_column+1):ncol(cover_data)]
  }


  total <- rowSums(data)

  shannon_index <- rowSums(-data[,1:ncol(data)]/total *log(data[,1:ncol(data)]/total), na.rm = TRUE)
  return(shannon_index)
}

hill_shannon <- function(cover_data, remove_column = NULL){
  if(is.null(remove_column)){

    vector <- shannon(cover_data)
    return(exp(vector))

  }
  else{
    vector <- shannon(cover_data, remove_column = remove_column)
    return(exp(vector))
  }


}

simpson <- function(cover_data, remove_column = NULL){

  if(is.null(remove_column)){
    data <- cover_data[,1:ncol(cover_data)]
  }
  else{
    n <- remove_column+1

    data <- cover_data[,n:ncol(cover_data)]
  }


  total <- rowSums(data)

  simpson_index <- rowSums((data[,1:ncol(data)]/total)^2)

  return(simpson_index)
}

hill_simpson <- function(cover_data, remove_column = NULL){
  if(is.null(remove_column)){

    vector <- simpson(cover_data)
    return(1/vector)

  }
  else{
    vector <- simpson(cover_data, remove_column = remove_column)
    return(1/vector)
  }


}


#A function to plot the different ways of measuring diversity against a gradient, i.e ph value,  for each plot.
plot_diversity <- function(data, diversities, plot_info, description = NULL) {
  plot_data <-  gather(data, key = "type", value = "Diversity", diversities)
  ggplot(data = plot_data, mapping = aes(x = as.numeric(plot_data[[plot_info]]), y = Diversity)) +
    geom_point()+
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(vars(type),scales = "free_y") +
    xlab(description) +
    ylab("Diversity estimates") +
    theme_update(strip.placement = "outside", strip.text.x = element_text(face = "bold"), title = element_text(face = "bold"))
}





beta <- function(cover, freq) {
  library(fitdistrplus)
  cover_data <- cover[,4:ncol(cover)]
  freq_data <- freq[,4:ncol(freq)]


  for (ele in colnames(cover_data)){
    if (sum(freq_data[,ele]) == 0){
      cover = cover[,!(names(cover) %in% ele)]
      cover_data <- cover_data[,!(names(cover_data) %in% ele)]
    }

  }

  #create data frame to hold the fitted values for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))

  colnames(beta_fit) <- c("species","a", "b")

  # for Each species calculate the shape parameter for the fitted beta distribution and save them in a data frame.
  for (specie in colnames(cover_data)) {
    beta_data <- cover_data[,specie]/16

    #remove all plots with 0 in frekvens.
    beta_data <- beta_data[freq_data[[specie]] == 1]


    if (length(unique(beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
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

  ggplot(data = df, aes(x = l)) +
    geom_line(aes(y = Observed, colour = "Observed cover data"))+
    geom_line(aes(y = Updated, colour = "New dataset"))+
    scale_colour_manual("",
                        values = c("Observed cover data"="red",
                                   "New dataset"="blue")) +

    labs(y = "Diversity", x = "Exponent l in Hill diveristy formula")+
    ggtitle(sprintf("Comparison of diversity estimates for plot %d", plot))




}

different_diversities2 <- function(data_observed, data_new, remove_column = NULL){
  observed <- data_observed[,(remove_column+1):ncol(data_observed)]
  update <- data_new[,(remove_column+1):ncol(data_new)]

  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("l", "Observed", "Updated")

  for (plot in 1:nrow(data_observed)) {
    for (l in seq(-1, 1, 0.1)){
    df[nrow(df)+1, ] <- c(l, diversity_index(observed[plot,],l = l) , diversity_index(update[plot,], l = l) )
    }
  }
  # Aggregating for all the different point
  df <- aggregate(df, list(df$l), mean )
  points <- df[df$l %in% c(-1,0,1),]
  ggplot(data = df, aes(x = l)) +
    geom_line(aes(y = Observed, colour = "Observed cover data"))+
    geom_line(aes(y = Updated, colour = "New dataset"))+
    geom_point(data = points, mapping = aes(x = as.numeric(l), y = Updated),
               fill = "blue", shape=15,  size = 2, colour = "blue") +
    geom_point(data = points, mapping = aes(x = as.numeric(l), y = Observed),
               fill = "red", shape=17, size = 2, colour = "red") +
    scale_colour_manual("",
                        values = c("Observed cover data"="red",
                                   "New dataset"="blue")) +

    labs(y = "Diversity", x = "Exponent l in Hill diveristy formula")+
    ggtitle(sprintf("Comparison of diversity estimates for %d plots", plot))


}
