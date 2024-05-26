
######################################################################################################
######################################## Sigma matrix ################################################
######################################################################################################

par(mfrow = c(4,4))
hist(sigma11, main = "", col = "lightgrey", xlab = "", ylab = "", xlim = range(0,15), probability = TRUE, breaks = 50,yaxt='n')
hist(corr12, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr13, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr14, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr12, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(sigma22, main = "", col = "lightgrey", xlab = "", ylab = "", xlim = range(0,15), probability = TRUE, breaks = 50,yaxt='n')
hist(corr23, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr24, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr13, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr23, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(sigma33, main = "", col = "lightgrey", xlab = "", ylab = "", xlim = range(0,15), probability = TRUE, breaks = 50,yaxt='n')
hist(corr34, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr14, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr24, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(corr34, main = "", xlim = range(-1,1),col = "lightgrey", xlab = "", ylab = "", probability = TRUE, breaks = 50,yaxt='n')
hist(sigma44, main = "", col = "lightgrey", xlab = "", ylab = "", xlim = range(0,15), probability = TRUE, breaks = 50,yaxt='n')



######################################################################################################
########################################### MCMC #####################################################
######################################################################################################

alpha_n <- length(posteriors$mu.alpha[1, ])
alpha_description <- c("edges", "gwesp", "gwdegree", "gwdsp")
mu.alpha <- matrix(rep(0, alpha_n),alpha_n,1)
par(mfrow=c(5,3))
for (k in 1:alpha_n) {
  # Constructing the dynamic y-label expression
  ylab_expr <- bquote(alpha[.(k)](.(alpha_description[k])))  # This adds description from alpha_description
  
  # Plotting the density plot for the current k
  alpha_vector <- as.vector(t(posteriors$mu.alpha[, k]))
  # Plot density
  plot(density(alpha_vector[seq(50, row_num, by = 50)]), lwd=2, main = "",
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
# Eta MCMC plots
Eta <- matrix(rep(0, 1),nrow = 1,1) #-1
## Plot density
eta_vector <- as.vector(t(posteriors$Eta))
plot(density(eta_vector[seq(50, row_num, by = 50)]), lwd=2,main = "",
     xlab=expression(eta), ylab="Density")
## Iterations
plot(ts(posteriors$Eta),ylim=range(c(posteriors$Eta,Eta)), main = "",
     ylab=expression(eta),xlab='iteration')
lines( c(1, dim(posteriors$Alpha)[1]),c(Eta,Eta), col='blue')
## Plot autocorrilations
acf(eta_vector, lag.max=1000, main = "",
    ylab="Autocorrelation", xlab="Lag")

######################################################################################################
############################### Alpha Interactions ###################################################
######################################################################################################



regular_length <- length(edges)
mu_length <- length(mu_edges)
# Create dataframe for plotting
df <- data.frame(edges = c(edges, mu_edges),
                 gwesp = c(gwesp, mu_gwesp),
                 gwdeg = c(gwdeg, mu_gwdeg),
                 group = c(rep("alpha", regular_length), rep("mu_alpha", mu_length)))

p1 <- ggplot(df, aes(x = edges, y = gwesp)) + 
  geom_point(aes(shape = group, color = group), size = 3) +
  scale_shape_manual(values = c("alpha" = 1, "mu_alpha" = 4)) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("alpha" = "black", "mu_alpha" = "red")) +
  #labs(title = "Plot of GWESP vs GWDEG", x = "GWESP", y = "GWDEG") +
  theme_classic()

p2 <- ggplot(df, aes(x = edges, y = gwdeg)) + 
  geom_point(aes(shape = group, color = group), size = 3) +
  scale_shape_manual(values = c("alpha" = 1, "mu_alpha" = 4)) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("alpha" = "black", "mu_alpha" = "red")) +
  #labs(title = "Plot of GWESP vs GWDEG", x = "GWESP", y = "GWDEG") +
  theme_classic()

p3 <- ggplot(df, aes(x = gwesp, y = gwdeg)) + 
  geom_point(aes(shape = group, color = group), size = 3) +
  scale_shape_manual(values = c("alpha" = 1, "mu_alpha" = 4)) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("alpha" = "black", "mu_alpha" = "red")) +
  #labs(title = "Plot of GWESP vs GWDEG", x = "GWESP", y = "GWDEG") +
  theme_classic()

ggarrange(p1, p2, p3,
          ncol = 3, nrow = 1,common.legend = TRUE, legend = "bottom")

######################################################################################################
######################### Alpha Granular Interaction Plots ###########################################
######################################################################################################

names <-  c("Siren (C)", "Cielnet (C)", "Acero (C)", "Jake (C)",
            "Juanes (C)", "Mambo (C)", "Natarajan (C)", "Togo (C)",
            "WT1 (C)", "WT2 (C)", "WT3 (C)", "Mali (T)",
            "Noordin Top (T)", "Al Quaeda Post (T)", "Al Quaeda Pre (T)", 
            "IS Europe (T)", "Ergekon (T)")

alpha <- full_model$Alpha
# Select all mu_alphas except gwdsp
mu_alpha <- full_model$mu.alpha[, 1:3]
colnames(mu_alpha) <- c("edges", "gwesp", "gwdeg")
row_num <- nrow(alpha)
#Extract every 100th row for thinning
# edges <- alpha[, 1, ][seq(100, row_num, by = 100), ]
# gwesp <- alpha[, 2, ][seq(100, row_num, by = 100), ]
# gwdeg <- alpha[, 3, ][seq(100, row_num, by = 100), ]
# mu_edges <- mu_alpha[,1][seq(100, row_num, by = 100)]
# mu_gwesp <- mu_alpha[,2][seq(100, row_num, by = 100)]
# mu_gwdeg <- mu_alpha[,3][seq(100, row_num, by = 100)]

#Extract every 300th row for thinning
edges <- alpha[, 1, ][seq(300, row_num, by = 300), ]
gwesp <- alpha[, 2, ][seq(300, row_num, by = 300), ]
gwdeg <- alpha[, 3, ][seq(300, row_num, by = 300), ]
mu_edges <- mu_alpha[,1][seq(300, row_num, by = 300)]
mu_gwesp <- mu_alpha[,2][seq(300, row_num, by = 300)]
mu_gwdeg <- mu_alpha[,3][seq(300, row_num, by = 300)]

# setting up for binding with names
colnames(edges) <- names
colnames(gwesp) <- names
colnames(gwdeg) <- names

edges <- as.data.frame(cbind(edges, mu_edges)) %>% 
  pivot_longer(everything(), names_to = "Group", values_to = "Edges")
gwesp <- as.data.frame(cbind(gwesp, mu_gwesp)) %>% 
  pivot_longer(everything(), names_to = "Group", values_to = "Gwesp")
gwdeg <- as.data.frame(cbind(gwdeg, mu_gwdeg)) %>% 
  pivot_longer(everything(), names_to = "Group", values_to = "Gwdeg")

df <- cbind(edges, gwesp$Gwesp, gwdeg$Gwdeg) 
colnames(df) <- c("Group", "Edges", "Gwesp", "Gwdeg")

other_categorisation <- c("Siren (C)", "Cielnet (C)", "Acero (C)",
                          "Juanes (C)", "Mambo (C)", "Natarajan (C)", "Togo (C)",
                          "WT3 (C)", "Mali (T)",
                          "Noordin Top (T)")


df$Group[df$Group  %in%  other_categorisation] <- "Other"
df$Group[df$Group=="mu_edges"] <- "mu Alpha"

names <-  c("Siren (C)", "Cielnet (C)", "Acero (C)", "Jake (C)",
            "Juanes (C)", "Mambo (C)", "Natarajan (C)", "Togo (C)",
            "WT1 (C)", "WT2 (C)", "WT3 (C)", "Mali (T)",
            "Noordin Top (T)", "Al Quaeda Post (T)", "Al Quaeda Pre (T)", 
            "IS Europe (T)", "Ergekon (T)")


p1 <- ggplot(df, aes(x = Edges, y = Gwesp)) + 
  geom_point(aes(shape = Group, color = Group), size = 3) +
  scale_shape_manual(values = c("mu Alpha" = 4, "Other" = 1, "Jake (C)" = 1, "WT1 (C)" = 1, "WT2 (C)" = 1,
                                "Al Quaeda Post (T)" = 1,  "Al Quaeda Pre (T)" = 1,  "IS Europe (T)" = 1,
                                "Ergekon (T)" = 1),
  ) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("mu Alpha" = "red", "Other" = "black", "Jake (C)" = "turquoise", "WT1 (C)" = "blue",
                                "WT2 (C)" = "blue", "Al Quaeda Post (T)" = "green",  "Al Quaeda Pre (T)" = "green",
                                "IS Europe (T)" = "orange", "Ergekon (T)" = "magenta")) +
  labs(x = "edges", y = "gwesp") +
  theme_classic()
p2 <- ggplot(df, aes(x = Edges, y = Gwdeg)) + 
  geom_point(aes(shape = Group, color = Group), size = 3) +
  scale_shape_manual(values = c("mu Alpha" = 4, "Other" = 1, "Jake (C)" = 1, "WT1 (C)" = 1, "WT2 (C)" = 1,
                                "Al Quaeda Post (T)" = 1,  "Al Quaeda Pre (T)" = 1,  "IS Europe (T)" = 1,
                                "Ergekon (T)" = 1),
  ) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("mu Alpha" = "red", "Other" = "black", "Jake (C)" = "turquoise", "WT1 (C)" = "blue",
                                "WT2 (C)" = "blue", "Al Quaeda Post (T)" = "green",  "Al Quaeda Pre (T)" = "green",
                                "IS Europe (T)" = "orange", "Ergekon (T)" = "magenta")) +
  labs( x = "edges", y = "gwdeg") +
  theme_classic()
p3 <- ggplot(df, aes(x = Gwesp, y = Gwdeg)) + 
  geom_point(aes(shape = Group, color = Group), size = 3) +
  scale_shape_manual(values = c("mu Alpha" = 4, "Other" = 1, "Jake (C)" = 1, "WT1 (C)" = 1, "WT2 (C)" = 1,
                                "Al Quaeda Post (T)" = 1,  "Al Quaeda Pre (T)" = 1,  "IS Europe (T)" = 1,
                                "Ergekon (T)" = 1),
  ) + # 1 is an empty circle, 4 is a cross
  scale_color_manual(values = c("mu Alpha" = "red", "Other" = "black", "Jake (C)" = "turquoise", "WT1 (C)" = "blue",
                                "WT2 (C)" = "blue", "Al Quaeda Post (T)" = "green",  "Al Quaeda Pre (T)" = "green",
                                "IS Europe (T)" = "orange", "Ergekon (T)" = "magenta")) +
  labs(x = "gwesp", y = "gwdeg") +
  theme_classic()


ggarrange(p1, p2, p3,
          ncol = 3, nrow = 1,common.legend = TRUE, legend = "bottom")

######################################################################################################
############################### Catterpillar plots ###################################################
######################################################################################################
par(mfrow = c(2,2))

parameters <- c("edges", "gwesp", "gwdegree", "gwdsp") #"log(n)"
k <- 1
#par(mfrow = c(2,2), las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

par(las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

network_order <- network_order %>% 
  mutate(Type = if_else(Number < 12, "Criminal", "Terrorist"))

mu_alpha <- posteriors$mu.alpha
for (k in c(1:4)){
  alpha <- posteriors$Alpha[,k,]
  colnames(alpha) <- network_order$Name
  # Calculating the medians for each column
  medians <- apply(alpha, 2, median)
  # Ordering the columns by median values
  order_indices <- order(medians)
  # Reordering the matrix columns based on medians and 
  alpha_ordered <- alpha[, order_indices]
  
  # Reorder according to medians
  network_order_for_col <- network_order[order(medians[network_order$Name]), ]
  
  # Assign colors based on 'Type'
  #colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
  # Creating the boxplot with reordered columns
  boxplot(alpha_ordered, 
          # col = colors,           # Use colors based on 'Type'
          col = "white",
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,
          ylab = paste(parameters[k]),
          #border = "black",
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          srt = 90)
  mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
  abline(h = mu_alpha99[1], lty = "dotted")
  abline(h = mu_alpha99[2], lty = "dotted")
  mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
  abline(h = mu_alpha90[1], lty = "dashed")
  abline(h = mu_alpha90[2], lty = "dashed")
}

######################################################################################################
########################### Catterpillar plots with N ################################################
######################################################################################################
par(mfrow = c(2,2))
networks <- data.frame(Number = 1:17, 
                       temp = c("Siren (C", "Cielnet (C", "Acero (C", "Jake (C",
                                "Juanes (C", "Mambo (C", "Natarajan (C", "Togo (C",
                                "WT1 (C", "WT2 (C", "WT3 (C", "Mali (T",
                                "Noordin Top (T", "Al Quaeda Post (T", "Al Quaeda Pre (T", 
                                "IS Europe (T", "Ergekon (T"),
                       N = unlist(lapply(nets, network.size)))
networks <- networks %>% 
  mutate(Name = paste(temp,",",  N, ")", sep = "")) %>% 
  select(Number, Name)


parameters <- c("edges", "gwesp", "gwdegree", "gwdsp") #"log(n)"
#par(mfrow = c(2,2), las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

par(las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

mu_alpha <- posteriors$mu.alpha
for (k in c(1:4)){
  alpha <- posteriors$Alpha[,k,]
  colnames(alpha) <- networks$Name
  # Calculating the medians for each column
  medians <- apply(alpha, 2, median)
  # Ordering the columns by median values
  order_indices <- order(medians)
  # Reordering the matrix columns based on medians and 
  alpha_ordered <- alpha[, order_indices]
  
  # Reorder according to medians
  network_order_for_col <- networks[order(medians[networks$Name]), ]
  
  # Assign colors based on 'Type'
  #colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
  # Creating the boxplot with reordered columns
  boxplot(alpha_ordered, 
          # col = colors,           # Use colors based on 'Type'
          col = "white",
          outline = TRUE,          # Ensures that outliers are plotted
          outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
          outbg = NA,
          ylab = paste(parameters[k]),
          #border = "black",
          boxwex = 0.6,            # Adjust the width of the box plots
          whisklty = 2,            # Set whisker line type to dashed
          staplelty = 2, 
          srt = 90)
  mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
  abline(h = mu_alpha99[1], lty = "dotted")
  abline(h = mu_alpha99[2], lty = "dotted")
  mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
  abline(h = mu_alpha90[1], lty = "dashed")
  abline(h = mu_alpha90[2], lty = "dashed")
}

######################################################################################################
########################### Catterpillar plots with full vs H1 #######################################
######################################################################################################

k=1
parameters <- c("edges", "gwesp", "gwdegree", "gwdsp") #"log(n)"
#par(mfrow = c(2,2), las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

par(las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
par(mfrow = c(2,1))
# Full model
mu_alpha <- full_model$mu.alpha
alpha <- full_model$Alpha[,k,]
colnames(alpha) <- networks$Name
# Calculating the medians for each column
medians <- apply(alpha, 2, median)
# Ordering the columns by median values
order_indices <- order(medians)
# Reordering the matrix columns based on medians and 
alpha_ordered <- alpha[, order_indices]

# Reorder according to medians
network_order_for_col <- networks[order(medians[networks$Name]), ]

# Assign colors based on 'Type'
#colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
# Creating the boxplot with reordered columns
boxplot(alpha_ordered, 
        # col = colors,           # Use colors based on 'Type'
        col = "white",
        outline = TRUE,          # Ensures that outliers are plotted
        outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
        outbg = NA,
        ylab = paste(parameters[k]),
        #border = "black",
        boxwex = 0.6,            # Adjust the width of the box plots
        whisklty = 2,            # Set whisker line type to dashed
        staplelty = 2, 
        srt = 90)
mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
abline(h = mu_alpha99[1], lty = "dotted")
abline(h = mu_alpha99[2], lty = "dotted")
mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
abline(h = mu_alpha90[1], lty = "dashed")
abline(h = mu_alpha90[2], lty = "dashed")

# H3 model
mu_alpha <- h1_edge_model$mu.alpha
alpha <- h1_edge_model$Alpha[,k,]
colnames(alpha) <- networks$Name
# Calculating the medians for each column
medians <- apply(alpha, 2, median)
# Ordering the columns by median values
order_indices <- order(medians)
# Reordering the matrix columns based on medians and 
alpha_ordered <- alpha[, order_indices]

# Reorder according to medians
network_order_for_col <- networks[order(medians[networks$Name]), ]

# Assign colors based on 'Type'
#colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
# Creating the boxplot with reordered columns
boxplot(alpha_ordered, 
        # col = colors,           # Use colors based on 'Type'
        col = "white",
        outline = TRUE,          # Ensures that outliers are plotted
        outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
        outbg = NA,
        ylab = paste(parameters[k]),
        #border = "black",
        boxwex = 0.6,            # Adjust the width of the box plots
        whisklty = 2,            # Set whisker line type to dashed
        staplelty = 2, 
        srt = 90)
mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
abline(h = mu_alpha99[1], lty = "dotted")
abline(h = mu_alpha99[2], lty = "dotted")
mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
abline(h = mu_alpha90[1], lty = "dashed")
abline(h = mu_alpha90[2], lty = "dashed")



######################################################################################################
########################### Catterpillar plots with full vs H3 #######################################
######################################################################################################

k=2
parameters <- c("edges", "gwesp", "gwdegree", "gwdsp") #"log(n)"
#par(mfrow = c(2,2), las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))

par(las = 2, mar = c(7, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
par(mfrow = c(2,1))
# Full model
mu_alpha <- full_model$mu.alpha
alpha <- full_model$Alpha[,k,]
colnames(alpha) <- networks$Name
# Calculating the medians for each column
medians <- apply(alpha, 2, median)
# Ordering the columns by median values
order_indices <- order(medians)
# Reordering the matrix columns based on medians and 
alpha_ordered <- alpha[, order_indices]

# Reorder according to medians
network_order_for_col <- networks[order(medians[networks$Name]), ]

# Assign colors based on 'Type'
#colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
# Creating the boxplot with reordered columns
boxplot(alpha_ordered, 
        # col = colors,           # Use colors based on 'Type'
        col = "white",
        outline = TRUE,          # Ensures that outliers are plotted
        outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
        outbg = NA,
        ylab = paste(parameters[k]),
        #border = "black",
        boxwex = 0.6,            # Adjust the width of the box plots
        whisklty = 2,            # Set whisker line type to dashed
        staplelty = 2, 
        srt = 90)
mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
abline(h = mu_alpha99[1], lty = "dotted")
abline(h = mu_alpha99[2], lty = "dotted")
mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
abline(h = mu_alpha90[1], lty = "dashed")
abline(h = mu_alpha90[2], lty = "dashed")

# H3 model
mu_alpha <- h3_gwesp_model$mu.alpha
alpha <- h3_gwesp_model$Alpha[,k,]
colnames(alpha) <- networks$Name
# Calculating the medians for each column
medians <- apply(alpha, 2, median)
# Ordering the columns by median values
order_indices <- order(medians)
# Reordering the matrix columns based on medians and 
alpha_ordered <- alpha[, order_indices]

# Reorder according to medians
network_order_for_col <- networks[order(medians[networks$Name]), ]

# Assign colors based on 'Type'
#colors <- ifelse(network_order_for_col$Type == "Criminal", "white", "darkgrey")
# Creating the boxplot with reordered columns
boxplot(alpha_ordered, 
        # col = colors,           # Use colors based on 'Type'
        col = "white",
        outline = TRUE,          # Ensures that outliers are plotted
        outpch = 21,             # Type 21 is a filled circle; set bg to NA for unfilled
        outbg = NA,
        ylab = paste(parameters[k]),
        #border = "black",
        boxwex = 0.6,            # Adjust the width of the box plots
        whisklty = 2,            # Set whisker line type to dashed
        staplelty = 2, 
        srt = 90)
mu_alpha99 <- calculate_99_percentiles(mu_alpha[,k])
abline(h = mu_alpha99[1], lty = "dotted")
abline(h = mu_alpha99[2], lty = "dotted")
mu_alpha90 <- calculate_90_percentiles(mu_alpha[,k])
abline(h = mu_alpha90[1], lty = "dashed")
abline(h = mu_alpha90[2], lty = "dashed")


######################################################################################################
########################### GOF plots #######################################
######################################################################################################
sim_n <- 100
gofnets <- gof.bayes.hergm(form = form, nets = nets, posteriors=posteriors,MCMC.burnin=1000,nsim=sim_n)

# Simulation vs observed
par(mfrow = c(2,2))
for (k in 1:Num.Nets){
  degree_dist_gof(k=k, nets=nets, gofnets=gofnets, sim_n=sim_n)
  geodesic_dist_gof(k=k, nets=nets, gofnets=gofnets, sim_n=sim_n)
  edge_wise_partner_dist_gof(k=k, nets=nets, gofnets=gofnets, sim_n=sim_n)
  triad_cencus_dist_gof(k=k, nets=nets, gofnets=gofnets, sim_n=sim_n)
  mtext(paste(network_order[k, 2],": Bayesian GOF Analysis", sep = ""), side = 3, line = - 2, outer = TRUE)
}



