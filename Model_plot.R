library(dplyr)
library(ggplot2)
baad<-read.csv('baad_data.csv')

baad$group1 <- ifelse(grepl("A", baad$pft),'Angiosperm' ,
                      ifelse(grepl("G", baad$pft), "Gymnosperm", "NA"))

baadnew <- baad %>% filter(!is.na(pft))

baadnew$DH_m2<-(baadnew$d.bh)*baadnew$h.t
baadnew$D2_3_m<-(baadnew$d.bh)^(2/3)

########################################
profit<-stack("Net_carbon_profit_mean0.05_LAImax_growing.nc")##mol C m-2 a-1
ncpKg<-(profit*12.01)/1000
baadpoint <- cbind(baadnew$longitude, baadnew$latitude)

baadnew$Pn<-raster::extract(ncpKg,baadpoint,method = 'bilinear')

baadnew <- baadnew %>%
  mutate(predicted_D = case_when(
    group1 == "Angiosperm" ~ (3 * 8.31 * Pn) / (8 * 177.74),
    group1 == "Gymnosperm" ~ (3 * 1.82 * Pn) / (8 * 92.09)
  ))


baadnew$dD<-baadnew$d.bh/baadnew$age
Ddata<-data.frame(baadnew[,c(3,4,17,36,63,67,68)])
Ddata <- Ddata %>%
  filter(!is.na(d.bh) & !is.na(age))

process_group <- function(df) {
  df %>%
    group_by(longitude, latitude, group1) %>%  # Group by location and plant group
    summarise(
      total_dD_dt = sum(dD, na.rm = TRUE),   # Calculate total dD/dt for each group
      total_predicted_D = sum(predicted_D, na.rm = TRUE),  # Sum predicted_D for each group
      .groups = 'drop'
    ) %>%
    filter(total_dD_dt != 0) %>%  # Exclude rows where total_dD_dt is 0%>%
    filter(total_predicted_D != 0) 
}

processed_data <- process_group(Ddata)

# Process data for Angiosperm group
processed_angiosperm <-  processed_data%>%
  filter(group1 == "Angiosperm") 

# Calculate median values for Gymnosperm group
processed_gymnosperm <-  processed_data%>%
  filter(group1 == "Gymnosperm") 


Dtotal_group1lm_inter <- lm(log(total_dD_dt) ~ log(total_predicted_D), data = processed_angiosperm)
summary(Dtotal_group1lm_inter)

Dtotal_group2lm_inter <- lm(log(total_dD_dt) ~ log(total_predicted_D), data = processed_gymnosperm)
summary(Dtotal_group2lm_inter)


group_colors <- c("Angiosperm" = "dodgerblue3", "Gymnosperm" = "darkorange")

lntotalDplotlm<-ggplot(plotdeta,aes(log(total_predicted_D), log(total_dD_dt),color = group1))+
  geom_point(size = 1) +
  xlab("lnModelled dD/dt")+
  ylab("lnObserved dD/dt")+
  xlim(-8.1,3)+ylim(-8.1,3)+
  #labs(color = "PFT") +
  #stat_quantile(quantiles = 0.90, formula = y ~ x - 1, size = 1, linetype = "solid") +
  #stat_quantile(quantiles = 0.90, formula = y ~ x - 1, size = 1, linetype = "dashed") +
  geom_smooth(aes(group = group1, color = group1),method = "lm",formula=y~x, na.rm=T,se = T) +
  geom_abline(slope=1, intercept = 0,col='forestgreen', linetype = "dashed",size = 1.5)+
  scale_color_manual(values = group_colors) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.text = element_text(size = 13),  # Adjust legend text size
    legend.title = element_blank()
  )
lntotalDplotlm

########################################################biomass
baadnew <- baadnew %>%
  mutate(predicted_biomass = case_when(
    group1 == "Angiosperm" ~ (177.74/8.31)*d.bh,
    group1 == "Gymnosperm" ~ (92.09/1.82)*d.bh
  ))


baadnew$observed_biomass<-(baadnew$m.to)*0.47/baadnew$a.cp

Bdata<-data.frame(baadnew[,c(3,4,31,50,63,69,70)])
Bdata <- Bdata %>%
  filter(!is.na(m.to) & !is.na(a.cp)& !is.na(predicted_biomass))


# Process data for Angiosperm group
processed_angiospermB <-  Bdata%>%
  filter(group1 == "Angiosperm") 

# Calculate median values for Gymnosperm group
processed_gymnospermB <-  Bdata%>%
  filter(group1 == "Gymnosperm") 



Btotal_group1_inter <- lm(log(observed_biomass) ~ log(predicted_biomass), data = processed_angiospermB)
summary(Btotal_group1_inter)##


Btotal_group2 <- lm(log(observed_biomass) ~ log(predicted_biomass) - 1, data = processed_gymnospermB)
summary(Btotal_group2)## 

Btotal_group2_inter <- lm(log(observed_biomass) ~ log(predicted_biomass), data = processed_gymnospermB)
summary(Btotal_group2_inter)## 

combined_plotBlm <- ggplot(Bdata, aes(log(predicted_biomass), log(observed_biomass), color = group1)) +
  geom_point(size = 1) +
  xlab(expression(lnModelled~Wtot~(KgC/m^2)))+
  ylab(expression(lnObserved~Wtot~(KgC/m^2)))+
  xlim(-2.1,4.1)+ylim(-2.1,4.1)+
  #stat_quantile(quantiles = 0.90, formula = y ~ x - 1, size = 1, linetype = "solid") +
  #stat_quantile(quantiles = 0.10, formula = y ~ x - 1, size = 1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = TRUE) +
  geom_abline(slope=1, intercept = 0,col='forestgreen', linetype = "dashed",size = 1)+
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 10)) +
  scale_color_manual(values = group_colors) +
  #labs(color = "PFT") +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )

combined_plotBlm

#######################################################
#######################################################turnover
baadnew <- baadnew %>%
  mutate(predict_turnover = case_when(
    group1 == "Angiosperm" ~ ((8*177.74)/(3*8.31))*((baadnew$d.bh)/(baadnew$Pn)),
    group1 == "Gymnosperm" ~ ((8*92.09)/(3*1.82))*((baadnew$d.bh)/(baadnew$Pn))
  ))

GPP<-stack("pythonGPP0.05_sm.nc") ##gc/m2/d
GPP <- calc(GPP, fun = function(x) mean(x, na.rm = TRUE))
plot(GPP)
NPP<-GPP*0.46*365/1000 ##KgC/m2/y
plot(NPP)
baadnew$NPP<-raster::extract(NPP,baadpoint,method = 'bilinear')

baadnew$observed_turnover<-baadnew$observed_biomass/baadnew$NPP


Tdata<-data.frame(baadnew[,c(3,4,63,71,73)])
Tdata <- Tdata %>%
  filter(!is.na(predict_turnover) & !is.na(observed_turnover))


# Process data for Angiosperm group
processed_angiospermT <-  Tdata%>%
  filter(group1 == "Angiosperm") 

# Calculate median values for Gymnosperm group
processed_gymnospermT <-  Tdata%>%
  filter(group1 == "Gymnosperm") 


Ttotal_group1_inter <- lm(log(observed_turnover) ~ log(predict_turnover), data = processed_angiospermT)
summary(Ttotal_group1_inter)##

Ttotal_group2_inter <- lm(log(observed_turnover) ~ log(predict_turnover), data = processed_gymnospermT)
summary(Ttotal_group2_inter)## 

combined_plotTlm <- ggplot(Tdata, aes(log(predict_turnover), log(observed_turnover), color = group1)) +
  geom_point(size = 1) +
  xlim(-1.7,4.5)+ylim(-1.7,4.5)+
  xlab(expression(lnModelled~tau~(yr)))+
  ylab(expression(lnObserved~tau~(yr)))+
  #stat_quantile(quantiles = 0.90, formula = y ~ x - 1, size = 1, linetype = "solid") +
  #stat_quantile(quantiles = 0.10, formula = y ~ x - 1, size = 1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, se = TRUE) +
  geom_abline(slope=1, intercept = 0,col='forestgreen', linetype = "dashed",size = 1)+
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 10)) +
  scale_color_manual(values = group_colors) +
  labs(color = "PFT") +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )

combined_plotTlm

###########################################################################
############################ 95th percentile of tree heights at each site
###########################################################################
baad2<-baadnew
baad2$anew<-73.64
equation_new <- function(Hm, H, anew, D) {
  H - Hm * (1 - exp(-anew * D / Hm))
}

# Define the find_Hm function
find_Hm_new <- function(H, anew, D) {
  # Provide an initial guess for Hm based on the range of data
  initial_guess <- mean(c(H, H * (1 - exp(-anew * D / H))))
  
  # Find the root using uniroot
  root <- try(uniroot(equation_new, interval = c(0.01, 200), H = H, anew = anew, D = D, extendInt = "yes", tol = 0.0001), silent = TRUE)
  
  # Check if uniroot was successful
  if(inherits(root, "try-error")) {
    # If uniroot failed, return NA
    return(NA)
  } else {
    # Extract the solution for Hm
    Hm_solution <- root$root
    # Check if Hm_solution is less than H
    if(Hm_solution < H) {
      return(NA)  # If Hm_solution is less than H, return NA
    } else {
      return(Hm_solution)  # Otherwise, return Hm_solution
    }
  }
}
# Example usage
H <- baad2$h.t
anew <- baad2$anew
D <- baad2$d.bh 

Hm_solutions <- mapply(find_Hm_new, H, anew, D)

baad2$Hmax <- Hm_solutions

baad2$Hmax[baad2$Hmax >160] <- NA
baad2$Hmax[baad2$Hmax <baad2$h.t] <- NA

Hdata<-data.frame(baad2[,c(3,4,63,66,75)])
Hdata <- Hdata %>%
  filter(!is.na(Hmax) & !is.na(Pn))##1972

library(quantreg)

# Fit quantile regression model at the 95th percentile (τ = 0.95)
library(dplyr)
process_groupH <- function(df) {
  df %>%
    group_by(longitude, latitude, group1) %>%  # Group by location and plant group
    summarise(
      Hmax_mean = mean(Hmax, na.rm = TRUE),    # Mean of Hmax
      Pn_mean = mean(Pn, na.rm = TRUE),        # Mean of Pn
      .groups = 'drop'
    ) %>%
    filter(Hmax_mean != 0)  # Exclude rows where Hmax is zero (if needed)
}

processed_dataH <- process_groupH(Hdata)

# Process data for Angiosperm group
processed_angiospermH <-  processed_dataH%>%
  filter(group1 == "Angiosperm") 

# Calculate median values for Gymnosperm group
processed_gymnospermH <-  processed_dataH%>%
  filter(group1 == "Gymnosperm") 

# Assuming processed_dataH is your result from process_groupH()
quant_model_95 <- rq(Hmax_mean ~ Pn_mean-1, tau = 0.95, data = processed_dataH)
summary(quant_model_95,se = "boot")


hmaxmodel_group1<-rq(Hmax_mean ~ Pn_mean-1, tau = 0.95,  data = processed_angiospermH)
summary(hmaxmodel_group1,se = "boot")##

hmaxmodel_group2<-rq(Hmax_mean ~ Pn_mean-1, tau = 0.95, data = processed_gymnospermH)
summary(hmaxmodel_group2,se = "boot")##

# Fit quantile regression model (95th percentile)
model <- rq(Hmax_mean ~ Pn_mean - 1, tau = 0.95, data = processed_dataH)

# Create new data for prediction
newdata <- data.frame(Pn_mean = seq(
  min(processed_dataH$Pn_mean, na.rm = TRUE),
  max(processed_dataH$Pn_mean, na.rm = TRUE),
  length.out = 100
))

# Predict with bootstrapped confidence intervals
boot_preds <- predict(model, newdata = newdata, se = "boot", interval = "confidence", level = 0.95)

# Combine predictions into data frame for plotting
plot_data <- cbind(newdata, as.data.frame(boot_preds))
colnames(plot_data) <- c("Pn_mean", "fit", "lwr", "upr")  # Rename for clarity

# Plot
ggplot() +
  # Scatter points from original data
  geom_point(data = processed_dataH, aes(x = Pn_mean, y = Hmax_mean, color = group1), size = 1) +
  
  # Quantile regression line and CI ribbon
  geom_line(data = plot_data, aes(x = Pn_mean, y = fit), color = "black", size = 1.2) +
  geom_ribbon(data = plot_data, aes(x = Pn_mean, ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.4) +
  
  # Labels and theme
  xlab(expression(Net~carbon~profit~(Kg~C/m^2~y))) +
  ylab("Hmax (m)") +
  scale_color_manual(values = group_colors) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )

############################add se for each group
# Set up prediction grid
newdata <- data.frame(Pn_mean = seq(
  min(processed_dataH$Pn_mean, na.rm = TRUE),
  max(processed_dataH$Pn_mean, na.rm = TRUE),
  length.out = 100
))

# ========== 1. Function to get bootstrapped CIs ==========
get_qr_ci <- function(model_data, newdata, tau = 0.95) {
  mod <- rq(Hmax_mean ~ Pn_mean - 1, tau = tau, data = model_data)
  pred <- predict(mod, newdata = newdata, se = "boot", interval = "confidence", level = 0.95)
  df <- cbind(newdata, as.data.frame(pred))
  colnames(df) <- c("Pn_mean", "fit", "lwr", "upr")
  return(df)
}

# ========== 2. Get CIs for overall, angiosperms, gymnosperms ==========
ci_all <- get_qr_ci(processed_dataH, newdata)
ci_angio <- get_qr_ci(filter(processed_dataH, group1 == "Angiosperm"), newdata)
ci_gymno <- get_qr_ci(filter(processed_dataH, group1 == "Gymnosperm"), newdata)

# Add group label for plotting
ci_angio$group <- "Angiosperm"
ci_gymno$group <- "Gymnosperm"

# Combine group-specific data
ci_groups <- bind_rows(ci_angio, ci_gymno)

# ========== 3. Final plot ==========
ggplot() +
  # Original data points
  geom_point(data = processed_dataH, aes(x = Pn_mean, y = Hmax_mean, color = group1), size = 1) +
  # Group-specific fits
  geom_line(data = ci_groups, aes(x = Pn_mean, y = fit, color = group), size = 1) +
  geom_ribbon(data = ci_groups, aes(x = Pn_mean, ymin = lwr, ymax = upr, fill = group), alpha = 0.2) +
  
  # Overall fit
  geom_line(data = ci_all, aes(x = Pn_mean, y = fit), color = "black", size = 1.2) +
  geom_ribbon(data = ci_all, aes(x = Pn_mean, ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.4) +
  
  
  # Labels and theme
  xlab(expression(Net~carbon~profit~(Kg~C/m^2~y))) +
  ylab("Hmax (m)") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )

cor(processed_angiospermH$Pn_mean, processed_angiospermH$Hmax_mean, use = "complete.obs", method = "pearson")

cor(processed_gymnospermH$Pn_mean, processed_gymnospermH$Hmax_mean, use = "complete.obs", method = "pearson")

cor(processed_dataH$Pn_mean, processed_dataH$Hmax_mean, use = "complete.obs", method = "pearson")

########################################## a global data set as the reference 
tallo<-read.csv("Tallo.csv")
tallo$D_m<-tallo$stem_diameter_cm/100
####################
talloheight<-data.frame(tallo[,cbind(6,7,9,14)])
talloheight<-na.omit(talloheight)
initial_Hmax <- max(talloheight$height_m, na.rm = TRUE)##115.8

lin_model <- lm(height_m ~ D_m-1, data = talloheight)##46.57548 
summary(lin_model)

initial_a <- coef(lin_model)[1] ##46.57548 

model <- nls(
  height_m ~ Hmax * (1 - exp(-a * D_m / Hmax)),
  data = talloheight,
  start = list(Hmax = initial_Hmax, a = initial_a),  # Initial guesses for parameters
  control = nls.lm.control(maxiter = 200) # Control parameters to handle iterations and warnings
)
summary(model)


# Extract the estimated coefficients
Hmax_est <- coef(model)["Hmax"]##33.2
a_est <- coef(model)["a"]##86.2476 

talloheight$predicted_height <- predict(model)
# Plot the observed data and the fitted curve
ggplot(talloheight, aes(x = D_m, y = height_m)) +
  geom_point(color = "gray", alpha = 0.8) +  # Plot the observed data points
  geom_line(aes(y = predicted_height), color = "black", size = 1)  +
  labs(
    x = "DBH (m)",
    y = "Height (m)"
  ) +
  theme_light() + theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+ annotate(
    "text", x = max(baadnewheight$d.bh) * 0.6, y = max(baadnewheight$h.t) * 0.1,
    label = paste0(
      "Hmax = ", round(Hmax_est, 2), "\n",
      "a = ", round(a_est, 2),"\n",
      "H = Hm(1 - exp(-aD/Hm))\n"
    ),
    hjust = 0
  )

###############################

# Define the objective function for Nelder-Mead optimization
objective_function <- function(a) {
  Hmax <- initial_Hmax  # Use the estimated or fixed Hmax value
  predicted_heights <- Hmax * (1 - exp(-a * talloheight$D_m / Hmax))  # Calculate predicted heights using the model
  sse <- sum((talloheight$height_m - predicted_heights)^2, na.rm = TRUE)    # Calculate the sum of squared errors
  return(sse)  # The goal is to minimize this value
}

# Initial guess for 'a'
initial_a <- initial_a # Set an initial guess close to what you've estimated before

# Optimize using Nelder-Mead method to find the best 'a'
opt_result <- optim(
  par = initial_a,                  # Initial value of 'a'
  fn = objective_function,          # The objective function defined above
  method = "Nelder-Mead",           # Derivative-free optimization method
  control = list(maxit = 1000)      # Control parameters, allowing sufficient iterations
)

# Extract the optimized 'a' value
optimized_a <- opt_result$par ##57.8

# Use optimized 'a' to refit the curve with fixed Hmax
Hmax <- initial_Hmax  # Assuming Hmax remains fixed for simplicity
fitted_heights <- Hmax * (1 - exp(-optimized_a * talloheight$D_m / Hmax))

# Plot the observed data and the fitted curve using the optimized 'a'
ggplot(talloheight, aes(x = D_m, y = height_m)) +
  geom_point(color = "gray", alpha = 0.8) +  # Observed data points
  geom_line(aes(y = fitted_heights), color = "black", size = 1) +  # Fitted curve with optimized 'a'
  labs(
    x = "DBH (m)",
    y = "Height (m)"
  ) +
  theme_light() + theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+
  annotate(
    "text", x = max(talloheight$D_m) * 0.6, y = max(talloheight$height_m) * 0.1,
    label = paste0(
      "Optimized parameters:\n",
      "Hmax = ", round(Hmax, 2), "\n",
      "a = ", round(optimized_a, 2),"\n",
      "H = Hm(1 - exp(-aD/Hm))\n"
    ),
    hjust = 0
  )

################
H <- tallo$height_m
tallo$anew<-57.76
anew <- tallo$anew
D <- tallo$D_m

Hm_solutionsTallo <- mapply(find_Hm_new, H, anew, D)
tallo$Hmax <- Hm_solutionsTallo

tallo$Hmax[tallo$Hmax >160] <- NA
tallo$Hmax[tallo$Hmax <tallo$height_m] <- NA

tallopoint <- cbind(tallo$longitude, tallo$latitude)
tallo$Pn<-raster::extract(ncpKg,tallopoint,method = 'bilinear')

process_tallo<- function(df) {
  df %>%
    group_by(longitude, latitude) %>%  # Group by location and plant group
    summarise(
      Hmax_mean = mean(Hmax, na.rm = TRUE),    # Mean of Hmax
      Pn_mean = mean(Pn, na.rm = TRUE),        # Mean of Pn
      .groups = 'drop'
    ) %>%
    filter(Hmax_mean != 0)  # Exclude rows where Hmax is zero (if needed)
}

site_tallo <- process_tallo(tallo)
Tallo_model_95 <- rq(Hmax_mean ~ Pn_mean-1, tau = 0.95, data = site_tallo)
summary(Tallo_model_95,se = "boot")

Tallo_model_95_test <- rq(Hmax_mean ~ Pn_mean, tau = 0.95, data = site_tallo)
summary(Tallo_model_95_test,se = "boot")
model <- rq(Hmax_mean ~ Pn_mean - 1, tau = 0.95, data = site_tallo)

# Create new data for prediction
newdata <- data.frame(Pn_mean = seq(
  min(site_tallo$Pn_mean, na.rm = TRUE),
  max(site_tallo$Pn_mean, na.rm = TRUE),
  length.out = 100
))

# Predict with bootstrapped confidence intervals
boot_preds <- predict(model, newdata = newdata, se = "boot", interval = "confidence", level = 0.95)

# Combine predictions into data frame for plotting
plot_data <- cbind(newdata, as.data.frame(boot_preds))
colnames(plot_data) <- c("Pn_mean", "fit", "lwr", "upr")  # Rename for clarity

# Plot
ggplot() +
  # Scatter points from original data
  geom_point(data = site_tallo, aes(x = Pn_mean, y = Hmax_mean), size = 1,alpha=0.5,color = "darkorchid") +
  
  # Quantile regression line and CI ribbon
  geom_ribbon(data = plot_data, aes(x = Pn_mean, ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
  geom_line(data = plot_data, aes(x = Pn_mean, y = fit), color = "darkorchid4", size = 1.2) +
  
  # Labels and theme
  xlab(expression(Net~carbon~profit~(Kg~C/m^2~y))) +
  ylab("Hmax (m)") +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_blank()
  )
cor(site_tallo$Pn_mean, site_tallo$Hmax_mean, use = "complete.obs", method = "pearson")
################################global Hm distribution
coefficient95 <- coef(Tallo_model_95)[1]##56.2
preHm95new<-ncpKg*coefficient95
plot(preHm95new)
