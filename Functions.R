# MCMC diagnotics function
mcmc_diagnostics <- function(posteriors){
  alpha_n <- length(posteriors$mu.alpha[1, ])
  alpha_description <- c("edges", "gwesp", "gwdegree", "gwdsp")
  mu.alpha <- matrix(rep(0, alpha_n),alpha_n,1)
  par(mfrow=c(1,3))
  for (k in 1:alpha_n) {
    # Constructing the dynamic y-label expression
    ylab_expr <- bquote(alpha[.(k)](.(alpha_description[k])))  # This adds description from alpha_description
    
    # Plotting the density plot for the current k
    alpha_vector <- as.vector(t(posteriors$mu.alpha[, k]))
    # Plot density
    plot(density(alpha_vector), lwd=2, main = "",
         xlab=ylab_expr, ylab="Density")
    
    # Plotting the graph
    plot(posteriors$Alpha[, k, 1], type='l', col='red',main = "",
         ylab=ylab_expr, ylim=range(posteriors$Alpha[, k, ]), xlab='iteration')
    
    # Adding lines for each network
    for (i in 1:Num.Nets) {
      lines(posteriors$Alpha[, k, i], col='red',)
    }
    
    # Additional lines as per original code
    lines(posteriors$mu.alpha[, k], col='black')
    lines(c(1, dim(posteriors$Alpha)[1]), c(mu.alpha[k], mu.alpha[k]), col='blue')
    
    # Plotting the autocorrelation plot for the current k
    acf(alpha_vector,  lag.max=1000, main = "",
        ylab="Autocorrelation", xlab="Lag")
  }
  #mtext("MCMC Output", side = 3, line = - 2, outer = TRUE)
  
  par(mfrow=c(1,3))
  # Eta MCMC plots
  Eta <- matrix(rep(0, 1),nrow = 1,1) #-1
  ## Plot density
  eta_vector <- as.vector(t(posteriors$Eta))
  plot(density(eta_vector), lwd=2,main = "",
       xlab=expression(eta), ylab="Density")
  ## Iterations
  plot(ts(posteriors$Eta),ylim=range(c(posteriors$Eta,Eta)), main = "",
       ylab=expression(eta),xlab='iteration')
  lines( c(1, dim(posteriors$Alpha)[1]),c(Eta,Eta), col='blue')
  ## Plot autocorrilations
  acf(eta_vector, lag.max=1000, main = "",
      ylab="Autocorrelation", xlab="Lag")
  
  
  # Beta MCMC plots
  Beta <- matrix(0,1,1) #-1
  ## Plot density
  beta_vector <- as.vector(t(posteriors$Beta[,1]))
  plot(density(beta_vector), lwd=2, main = "",
       xlab=expression(beta), ylab="Density")
  ## Iterations
  plot(ts(posteriors$Beta[,1]),ylim=range(c(posteriors$Beta[,1],Beta)), main = "",
       ylab=expression(beta),xlab='iteration')
  lines( c(1, dim(posteriors$Alpha)[1]),c(Beta,Beta), col='blue')
  ## Plot autocorrilations
  acf(beta_vector, lag.max=1000, main = "",
      ylab="Autocorrelation", xlab="Lag")
  
  #mtext("MCMC Output", side = 3, line = - 2, outer = TRUE)
}

# Define the function to calculate the desired percentiles
calculate_percentiles <- function(x) {
  quantile(x, probs = c(0.025, 0.975))
}

calculate_99_percentiles <- function(x) {
  quantile(x, probs = c(0.005, 0.995))
}
calculate_90_percentiles <- function(x) {
  quantile(x, probs = c(0.05, 0.95))
}

# Function for creating edge-wise shared partner dist of a network
shared_partners_fun <- function(net){
  summary(net ~ gwesp)
}

# Define the function to graph the simulated vs observed Degree Distribution
degree_dist_gof <- function(k, nets, gofnets, sim_n){
  obs_degree <- degree(nets[[k]], gmode = "graph")
  degrees <- degree(gofnets[[k]], g = c(1:sim_n), gmode = "graph")
  n <- network.size(nets[[k]])
  
  # Calculate the maximum degree including both observed and simulated networks
  max.deg <- max(c(obs_degree, degrees))
  
  # Create a matrix to hold the frequency of each degree from each simulation
  deg.sist <- cbind( matrix(colSums(degrees==0),sim_n,1),
                     t(apply(degrees,2,function(x) tabulate(x, nbins=max.deg) ) ) )# when tabulating we need to add the number of isolates 
  # Convert to densities
  deg.sist <- deg.sist/n
  # Convert matrix to data frame for easier handling with boxplot
  df <- as.data.frame(deg.sist)
  names(df) <- as.character(0:max.deg)
  
  # Obs degree freq table
  obs_freq <- as.data.frame(table(obs_degree))
  # Define the range of numbers
  range_of_numbers <- 0:(max.deg)
  # Create an empty dataframe to hold the final frequency table
  freq_table <- data.frame(number = range_of_numbers, frequency = 0)
  # Update the frequencies for existing numbers
  freq_table$frequency[freq_table$number %in% obs_freq$obs_degree] <- obs_freq$Freq
  
  # Fix observed densities
  density <- freq_table %>% 
    mutate(density = frequency/n)
  
  # Plotting boxplots with custom outlier symbols and custom 
  ylim <- c(0, max(density$density, deg.sist))
  boxplot(df, names = 0:max.deg,
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,              # Set background of outliers to NA for unfilled
          main = "",  #paste(network_order[k, 2]),
          xlab = "Degree",
          ylab = "",#"Density",
          border = "black",
          col = "white",       # Add some color to the boxplots
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          ylim = ylim)           # Set staple line type to dashed
  
  # Observed density distribution
  lines(density$number+1, density$density, 
        type = 'l', pch = 16, col = 'red', lwd = 2)
  
  # 97.5 and 2.5 percentile
  percentiles <- apply(deg.sist, 2, calculate_percentiles)
  
  lines(0:max.deg+1, percentiles[1,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
  lines(0:max.deg+1, percentiles[2,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
}

# Define the function to graph the simulated vs observed Shortest Geodesic Distance Distribution
geodesic_dist_gof <- function(k, nets, gofnets, sim_n){
  
  obs_geodist <- geodist(nets[[k]], count.paths = FALSE)$gdist
  
  # Use the table function to count the frequencies of each unique value
  obs_lower_triangle <- obs_geodist[lower.tri(obs_geodist, diag = FALSE)]
  obs_element_freq <- table(obs_lower_triangle)
  obs_density <- obs_element_freq/sum(obs_element_freq)
  
  # Creating the the density geodesic distance distributions from the sim graphs
  sim_graphs_list <- gofnets[[k]]
  sim_geodist <- lapply(sim_graphs_list, geodist, count.paths = FALSE)
  
  
  # Looping to extract the density geodesic distance distributions and the max
  sim_densities <- list()
  max_distances <- c()
  for (i in 1:sim_n){
    # Select the geodist matrix and calculate the distribution
    geodist_matrix <- sim_geodist[[i]]$gdist
    sim_lower_triangle <- geodist_matrix[lower.tri(geodist_matrix, diag = FALSE)]
    sim_element_freq <- table(sim_lower_triangle)
    sim_density <- sim_element_freq/sum(sim_element_freq)
    
    # Create list with all the geodesic distance distributions and vector with max distances
    sim_densities[[i]] <- sim_density
    
    string <- names(sim_density)
    string <- subset(string, string != "Inf")
    max_dist <- max(as.numeric(string))
    max_distances <- c(max_distances, max_dist)
  }
  
  #Creating the data frame with right columns
  max_dist <- max(max_distances)
  temp_matrix <- matrix(nrow = 1, ncol = max_dist+1)
  colnames(temp_matrix) <- c(as.character(1:max_dist), "Inf")
  df <- as.data.frame(temp_matrix)
  
  # Adding all the sim dist to the data frame
  for (i in 1:sim_n){
    temp <- as.data.frame(sim_densities[[i]]) %>% 
      pivot_wider(names_from = sim_lower_triangle, values_from = Freq)
    df <- df %>% 
      add_row(temp)
  }
  # Removing the first empty row, replacing all NAs with 0
  df <- df[-1,]
  df <- df %>% replace(is.na(.), 0) %>% 
    rename("NA" = "Inf")
  
  # Plotting boxplots with custom outlier symbols 
  ylim <- c(0, max(obs_density, max(df)))
  boxplot(df, 
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,              # Set background of outliers to NA for unfilled
          main =  "", #paste(network_order[k, 2]),
          xlab = "Min Geodesic Distance",
          ylab = "",#"Density",
          las = 2,                 # Make x-axis labels perpendicular for better readability
          border = "black",
          col = "white",       # Add some color to the boxplots
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          ylim = ylim, 
          # Remove axis as the lines cant handle non numeric arguments that would be due to NAs
          axes = FALSE)           
  
  # Observed distribution
  obs_df <- as.data.frame(temp_matrix)
  
  temp <- as.data.frame(obs_density) %>% 
    pivot_wider(names_from = obs_lower_triangle, values_from = Freq)
  obs_df <- obs_df %>% 
    add_row(temp)
  
  obs_df <- obs_df[-1,]
  obs_df <- obs_df %>% 
    replace(is.na(.), 0) %>% 
    rename("NA" = "Inf")
  
  lines(1:(max_dist+1), obs_df,
        type = 'l', pch = 16, col = 'red', lwd = 2)
  
  # 97.5 and 2.5 percentile
  percentiles <- sapply(df, calculate_percentiles)
  
  lines(1:(max_dist+1), percentiles[1,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
  lines(1:(max_dist+1), percentiles[2,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
  
  
  # Return the axis and border of graph
  label <- names(obs_df)
  # axis(1, 1:(max_dist+1))
  axis(2)
  axis(1, at = 1:length(label), labels = label)
  box(lwd = 1)
}

# Define the function to graph the simulated vs observed Triad Cencus
triad_cencus_dist_gof <- function(k, nets, gofnets, sim_n){
  
  # Generating triad census dist
  obs_triad <- triad.census(nets[[k]], mode = "graph")
  sim_triad <- triad.census(gofnets[[k]], g = 1:sim_n, mode = "graph")
  
  # converting to density
  total <- sum(obs_triad)
  obs_triad_df <- as.data.frame(obs_triad/total)
  sim_triad_df <- as.data.frame(sim_triad/total)
  
  
  # Plotting boxplots with custom outlier symbols and custom 
  ylim <- c(0, max(obs_triad_df, sim_triad_df))
  boxplot(sim_triad_df,
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,              # Set background of outliers to NA for unfilled
          main = "",  #paste(network_order[k, 2]),
          xlab = "Triad Cencus",
          ylab = "",#"Density",
          border = "black",
          col = "white",       # Add some color to the boxplots
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          ylim = ylim)           # Set staple line type to dashed
  
  # Observed density distribution
  lines(1:4, obs_triad_df, 
        type = 'l', pch = 16, col = 'red', lwd = 2)
  
  # 97.5 and 2.5 percentile
  percentiles <- apply(sim_triad/total, 2, calculate_percentiles)
  
  lines(1:4, percentiles[1,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
  lines(1:4, percentiles[2,], 
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
}

# Define the function to graph the simulated vs observed edge-wise shared partner
edge_wise_partner_dist_gof <- function(k, nets, gofnets, sim_n){
  # Calculate observed shared partner distribution, max obs and col number for temp matrix
  obs_shared_partners <- shared_partners_fun(nets[[k]])
  obs_max <- max(obs_shared_partners/sum(obs_shared_partners))
  if(is.na(obs_max)){
    obs_max <- 0
  }
  coln <- length(obs_shared_partners)
  
  # Calculate simulated shared partner distribution
  temp_matrix <- matrix(nrow = sim_n, ncol = coln)
  for (i in 1:sim_n){
    temp_matrix[i,] <- shared_partners_fun(gofnets[[k]][[i]])
  }
  
  # Find x-lim 
  dist <- colSums(temp_matrix)
  max <- 0
  for (i in coln:1){
    if(dist[i]!=0) {
      max <- i}
    if (max != 0) break
  }
  
  # Making it into densities
  vec <- apply(temp_matrix, 1, sum)
  density_matrix <- temp_matrix / vec
  # Fixing NaN
  density_matrix[is.na(density_matrix)] <- 0
  
  df <- as.data.frame(density_matrix)[,1:max]
  colnames(df) <- as.character(0:(max-1))
  
  
  
  # Plotting boxplots with custom outlier symbols and custom 
  ylim <- c(0, max(obs_max, density_matrix))
  boxplot(df,
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,              # Set background of outliers to NA for unfilled
          main = "", #paste(network_order[k, 2]),
          xlab = "Edge-Wise Shared Partners",
          ylab = "",#"Density",
          border = "black",
          col = "white",       # Add some color to the boxplots
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          ylim = ylim)           # Set staple line type to dashed
  
  # Observed density distribution
  obs_density <- (obs_shared_partners/sum(obs_shared_partners))[1:max]
  lines(1:max, obs_density,
        type = 'l', pch = 16, col = 'red', lwd = 2)
  
  # 97.5 and 2.5 percentile
  percentiles <- (apply(density_matrix, 2, calculate_percentiles))[,1:max]
  
  lines(1:max, percentiles[1,],
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
  lines(1:max, percentiles[2,],
        type = 'l', pch = 16, col = 'darkgray', lwd = 1)
}








