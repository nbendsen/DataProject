shanon_index_v2 <- function(cover, freq) {
  library(fitdistrplus)
  cover_data <- cover[,4:ncol(cover)]
  freq_data <- freq[,4:ncol(freq)]


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
  }

  # create list for shannon index for each plot
  shanon_list <- c()

  for (row in 1:nrow(cover_data)) {

    # Create an empty list for a given row
    mean_posterior <- c()


    # for a given row, find out what species is found in frekvens
    species_spotted_in_frekvens <- colnames(freq_data[c(freq_data[row,]  == 1)])

    #For each species spotten in frekvens, appends its posterior cover to the cover data for that row
    for (species_spotted in species_spotted_in_frekvens ) {


      alpha_post <- as.numeric((as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                                  as.numeric(cover_data[[species_spotted]][row]) ))
      beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + 16 - as.numeric(cover_data[[species_spotted]][row])

      mean_posterior <- append(mean_posterior, (alpha_post)/(alpha_post+beta_post))




    }

    #Calculate the shanon index value and append it to the list after normalizing and removing zeroes

    total_cover <- sum(mean_posterior)
    mean_posterior <- mean_posterior[!(mean_posterior < 0.00001)]
    shanon_value <- -sum(mean_posterior/total_cover * log((mean_posterior/total_cover)))

    shanon_list <- append(shanon_list,shanon_value)

  }
  out <- cover[,1:3]
  out$shannon <- shanon_list
  return(out)
}





#Function for Species richness from frequency data.
#the species indicator is a substring that is present in only the the columns
# with data for species. i.e not plot, site or year. If species indicator is left empty, all columns will be treated as species columns.

species_richness <- function(frequency_data ,species_indicator = NULL) {

  if (is.null(species_indicator)) {
    species_richness <- rowSums(frequency_data)
  }
  else {
    species_richness <- rowSums(frequency_data[, grep(species_indicator, names(frequency_data))])
  }
  return(species_richness)
}

#Function for the simple shannon index only using the cover data.
#the species indicator is a substring that is present in only the the columns
# with data for species. i.e not plot, site or year. If species indicator is left empty, all columns will be treated as species columns.

row_shannon <- function(cover_data ,species_indicator = NULL) {

  if (is.null(species_indicator)) {
    shannon <- rowSums(-cover_data/(rowSums(cover_data)) *log(cover_data/(rowSums(cover_data))), na.rm = TRUE  )
  }

  else {
    shannon <- rowSums(-cover_data[, grep(species_indicator, names(cover_data))]/
                         (rowSums(cover_data[,grep(species_indicator, names(cover_data))])) *
                         log(cover_data[, grep(species_indicator, names(cover_data))]/
                               (rowSums(cover_data[,grep(species_indicator, names(cover_data))]))), na.rm = TRUE  ) }
  return(shannon)
}


#A function to plot the different ways of measuring diversity against a gradient, i.e ph value,  for each plot.



plot_index_gradient <- function(data, indexes, gradient ,gradient_text = NULL) {
  long_data <-  gather(data, key = "type", value = "index", indexes)
  ggplot(data = long_data, mapping = aes(x = as.numeric(long_data[[gradient]]), y = index)) +
    geom_point()+
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(vars(type),scales = "free_y") +
    xlab(gradient_text) +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(face = "bold"), title = element_text(face = "bold"))
}




update_abundance_data <- function(abundance_data, present_data, n = 1 , remove_column = 0) {
  library(fitdistrplus)
  saved_columns <- abundance_data[,0:remove_column, drop = FALSE]
  abundance_data <- abundance_data[,(1 + remove_column):ncol(abundance_data)]
  present_data <- present_data[,(1 + remove_column):ncol(present_data)]


  #create data frame to hold the fitted values for each species
  beta_fit <- data.frame(matrix(ncol = 3, nrow = 0))

  colnames(beta_fit) <- c("species","a", "b")

  # for Each species calculate the shape parameter for the fitted beta distribution and save them in a data frame.
  for (specie in colnames(abundance_data)) {
    beta_data <- abundance_data[,specie]/n

    #remove all plots with 0 in frekvens.
    beta_data <- beta_data[present_data[[specie]] == 1]


    if (length(unique(beta_data)) > 1) {
      beta_data_fitted <- fitdist(beta_data, "beta", method = "mme")
      beta_fit[nrow(beta_fit) + 1,] <- c(specie, beta_data_fitted$estimate[1], beta_data_fitted$estimate[2])

    }
  }

  #output <-  data.frame(matrix(ncol = ncol(abundance_data), nrow = nrow(abundance_data)))
  #colnames(output) <- colnames(abundance_data[remove_column: ncol(abundance_data)])

  for (row in 1:nrow(present_data)) {

    # Create an empty list for a given row
    mean_posterior <- c()


    # for a given row, find out what species is found in present
    species_spotted_in_present <- colnames(present_data[c(present_data[row,]  == 1)])

    #For each species spotted in present, appends its posterior abundance_data to the abundance_data data for that row
    for (species_spotted in species_spotted_in_present ) {

      alpha_post <- as.numeric((as.numeric(beta_fit[beta_fit$species == species_spotted,]$a) +
                                  as.numeric(abundance_data[[species_spotted]][row]) ))
      beta_post <-  as.numeric(beta_fit[beta_fit$species == species_spotted,]$b) + n - as.numeric(abundance_data[[species_spotted]][row])

      if (length(alpha_post) == 1 && length(beta_post) == 1) {

        abundance_data[row, species_spotted] <- ((alpha_post)/(alpha_post+beta_post))
      }
      else if(!is.null(abundance_data[row, species_spotted]) && abundance_data[row, species_spotted] > 0) {
        abundance_data[row, species_spotted] <- 0
      }

    }
  }

  return(cbind(saved_columns ,abundance_data))
}



diversity_index <- function(abundance_data, remove_columns = 0, l = 1) {

  abundance_data <- abundance_data[,(1 + remove_columns):ncol(abundance_data)]


  shannon_list <- c()

  if (l == 1) {
    for (row in 1:nrow(abundance_data)) {

      all_obs <- abundance_data[row,]
      total_cover <- sum(all_obs)
      all_obs <- all_obs[(all_obs > 0.00001)]
      shannon_value <- -sum(all_obs/total_cover * log((all_obs/total_cover)))
      shannon_list <- append(shannon_list,shannon_value)
    }
  }
  else{
    for (row in 1:nrow(abundance_data)) {

      all_obs <- abundance_data[row,]
      total_cover <- sum(all_obs)
      all_obs <- all_obs[!(all_obs == 0)]
      shannon_value <- log((sum((all_obs/total_cover) * (1/(all_obs/total_cover))^l))^(1/l))
      shannon_list <- append(shannon_list, shannon_value)
    }
  }
  return(shannon_list)
}

