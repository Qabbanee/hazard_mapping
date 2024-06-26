library(INLA)
library(inlabru)
library(ggplot2)
library(readxl)
library(dplyr)

# Set bru options to compute DIC
bru_options_set(control.compute = list(dic = TRUE))

# Read avalanche data and filter by date
avalanche_data <- read_excel("~/Desktop/hazard_mapping/Lyngen_detections_attribute.xlsx")
avalanche_data$Dato <- as.Date(avalanche_data$Dato, format = "%Y-%m-%d")
start_date <- as.Date("2020-01-01")
end_date <- as.Date("2020-03-30")
avalanche_data <- avalanche_data[avalanche_data$Dato >= start_date & avalanche_data$Dato <= end_date, ]
avalanche_data <- avalanche_data[order(avalanche_data$skredTidspunkt), ]

# Read warning data
warning_data <- read_excel("~/Desktop/hazard_mapping/avalanche-bulletin_v1.xlsx")
warning_data$Dato <- as.Date(warning_data$Dato, format = "%Y-%m-%d")

# Prepare data for points with avalanches
point_dta <- data.frame(X = avalanche_data$X, Y = avalanche_data$Y, Date = avalanche_data$Dato)
point_dta$time <- as.integer(difftime(point_dta$Date, min(point_dta$Date), units = "days")) + 1

# Function to specify decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))

# Read and prepare elevation, slope, and aspect data
hoyde_data <- read_excel("~/Desktop/hazard_mapping/hoyde100m.xlsx")
slope_data <- read_excel("~/Desktop/hazard_mapping/slope100m.xlsx")
aspect_data <- read_excel("~/Desktop/hazard_mapping/aspect100m.xlsx")

elevData <- data.frame(X = hoyde_data$X, Y = hoyde_data$Y, elev = round(hoyde_data$LYNGEN100M_BAND_ / 1000, 2)) # in km
slopeData <- data.frame(X = slope_data$X, Y = slope_data$Y, slope = specify_decimal(slope_data$Slope_lyngen100m_Band_1, 1))
aspectData <- data.frame(X = aspect_data$X, Y = aspect_data$Y, aspect = specify_decimal(aspect_data$Aspect_lyngen100m_Band_1, 1))

# For factors, use all decimals
slopeData.Factor <- data.frame(X = slope_data$X, Y = slope_data$Y, slope = slope_data$Slope_lyngen100m_Band_1)
aspectData.Factor <- data.frame(X = aspect_data$X, Y = aspect_data$Y, aspect = aspect_data$Aspect_lyngen100m_Band_1)

# Convert data frames to spatial objects
coordinates(aspectData) <- ~X + Y
gridded(aspectData) <- TRUE
proj4string(aspectData) <- CRS("+proj=utm +zone=33 +datum=WGS84")

coordinates(slopeData) <- ~X + Y
gridded(slopeData) <- TRUE
proj4string(slopeData) <- CRS("+proj=utm +zone=33 +datum=WGS84")

coordinates(elevData) <- ~X + Y
gridded(elevData) <- TRUE
proj4string(elevData) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Functions to convert slope and aspect to factors
convert_to_compass_direction_2 <- function(value) {
  if (is.na(value)) return(NA)
  if (value == -1) return('A')
  if (value >= 0 && value < 22.5) return('N')
  if (value >= 22.5 && value < 67.5) return('NE')
  if (value >= 67.5 && value < 112.5) return('E')
  if (value >= 112.5 && value < 157.5) return('SE')
  if (value >= 157.5 && value < 202.5) return('S')
  if (value >= 202.5 && value < 247.5) return('SW')
  if (value >= 247.5 && value < 292.5) return('W')
  if (value >= 292.5 && value < 337.5) return('NW')
  return('N')
}

convert_grad <- function(value) {
  if (is.na(value)) return(NA)
  if (value >= 0 && value < 7) return('(0-7)°')
  if (value >= 7 && value < 13) return('(7-13)°')
  if (value >= 13 && value < 19) return('(13-19)°')
  if (value >= 19 && value < 25) return('(19-25)°')
  if (value >= 25 && value < 31) return('(25-31)°')
  if (value >= 31 && value < 37) return('(31-37)°')
  if (value >= 37 && value <= 45) return('(37-45)°')
  return('>45°')
}

aspectData.Factor$aspect <- sapply(aspectData.Factor$aspect, convert_to_compass_direction_2)
slopeData.Factor$slope <- sapply(slopeData.Factor$slope, convert_grad)
aspectData.Factor$aspect <- factor(aspectData.Factor$aspect)
slopeData.Factor$slope <- factor(slopeData.Factor$slope)

coordinates(aspectData.Factor) <- ~X + Y
gridded(aspectData.Factor) <- TRUE
proj4string(aspectData.Factor) <- CRS("+proj=utm +zone=33 +datum=WGS84")

coordinates(slopeData.Factor) <- ~X + Y
gridded(slopeData.Factor) <- TRUE
proj4string(slopeData.Factor) <- CRS("+proj=utm +zone=33 +datum=WGS84")

coordinates(point_dta) <- c("X", "Y")
proj4string(point_dta) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Define the boundary for the study area
utm_coords <- matrix(c(
  701695, 7778986,  # Point 1
  682052, 7740378,  # Point 2
  672061, 7737838,  # Point 3
  650857, 7723935,  # Point 4
  656483, 7707188,  # Point 5
  660716, 7700415,  # Point 6
  669816, 7699624,  # Point 7
  672134, 7695965,  # Point 8
  670426, 7688647,  # Point 9
  686526, 7683280,  # Point 10
  694453, 7691330,  # Point 11
  708235, 7718650,  # Point 12
  711406, 7738652,  # Point 13
  706762, 7772784,  # Point 14
  701695, 7778986
), ncol = 2, byrow = TRUE)

p <- Polygon(utm_coords)
ps <- Polygons(list(p), "boundary_id")
boundary <- SpatialPolygons(list(ps))
proj4string(boundary) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Create the mesh for spatial modeling
boundary_proj <- proj4string(boundary)
max.edge <- diff(range(coordinates(point_dta)[, 1])) / (3 * 5)
bound.outer <- diff(range(coordinates(point_dta)[, 1])) / 3

mesh <- inla.mesh.2d(
  loc = coordinates(point_dta),
  boundary = boundary,
  max.edge = c(1, 5) * max.edge,
  offset = c(max.edge, bound.outer),
  cutoff = max.edge / 3,
  crs = CRS(boundary_proj)
)

matern <- inla.spde2.pcmatern(
  mesh,
  prior.sigma = c(0.1, 0.01),
  prior.range = c(10, 0.1)
)

# Define and fit multiple models
cmp1S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "linear")
fit1S <- lgcp(cmp1S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp2S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T)))
fit2S <- lgcp(cmp2S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp3S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "linear")
fit3S <- lgcp(cmp3S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp4S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T)))
fit4S <- lgcp(cmp4S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp5S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(1000), fixed = T))) +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T))) +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast")
fit5S <- lgcp(cmp5S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp6S <- coordinates ~ Intercept(1) +
  mySmooth(coordinates, model = matern) +
  slope(f.slope(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(1000), fixed = T))) +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T))) +
  aspect(f.aspect(.data.), model = "rw1", cyclic = T, hyper = list(theta = list(initial = log(1000), fixed = T)))
fit6S <- lgcp(cmp6S, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp1 <- coordinates ~ Intercept(1) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "linear")
fit1 <- lgcp(cmp1, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp2 <- coordinates ~ Intercept(1) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T)))
fit2 <- lgcp(cmp2, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp3 <- coordinates ~ Intercept(1) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "linear")
fit3 <- lgcp(cmp3, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp4 <- coordinates ~ Intercept(1) +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T)))
fit4 <- lgcp(cmp4, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp5 <- coordinates ~ Intercept(1) +
  slope(f.slope(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(1000), fixed = T))) +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T))) +
  aspect(f.aspect.Factor(.data.), model = "factor_contrast")
fit5 <- lgcp(cmp5, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

cmp6 <- coordinates ~ Intercept(1) +
  slope(f.slope(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(1000), fixed = T))) +
  elev(f.elev(.data.), model = "rw1", cyclic = F, hyper = list(theta = list(initial = log(25), fixed = T))) +
  aspect(f.aspect(.data.), model = "rw1", cyclic = T, hyper = list(theta = list(initial = log(1000), fixed = T)))
fit6 <- lgcp(cmp6, data = point_dta, samplers = boundary, domain = list(coordinates = mesh))

# Compare models using WAIC and DIC criteria
knitr::kable(deltaIC(fit1S, fit2S, fit3S, fit4S, fit5S, fit6S, fit1, fit2, fit3, fit4, fit5, fit6, criterion = c("WAIC", "DIC")))

# Generate predictions for the models with mySmooth spatial effect
pred.df <- fm_pixels(mesh, mask = boundary, format = "sp")
fit.int1 <- predict(fit1, pred.df, ~ exp(mySmooth + elev + Intercept + aspect + slope))
fit.int2 <- predict(fit2, pred.df, ~ exp(mySmooth + elev + Intercept + aspect + slope))
fit.int3 <- predict(fit3, pred.df, ~ exp(mySmooth + elev + Intercept + slope))
fit.int4 <- predict(fit4, pred.df, ~ exp(mySmooth + elev + Intercept + slope))
fit.int5 <- predict(fit5, pred.df, ~ exp(mySmooth + elev + Intercept + aspect + slope))
fit.int6 <- predict(fit6, pred.df, ~ exp(mySmooth + elev + Intercept + aspect + slope))

# Plot the predictions
ggplot() +
  gg(fit.int1, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

ggplot() +
  gg(fit.int2, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

ggplot() +
  gg(fit.int3, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

ggplot() +
  gg(fit.int4, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

ggplot() +
  gg(fit.int5, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

ggplot() +
  gg(fit.int6, aes(fill = mean)) +
  gg(boundary, alpha = 0, lwd = 0.2) +
  scale_fill_continuous(name = "Mean Value", limits = c(0, max(fit.int1$mean))) +
  labs(x = "X", y = "Y") +
  theme_minimal()

