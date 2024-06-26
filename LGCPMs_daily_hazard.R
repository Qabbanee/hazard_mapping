library(INLA)
library(inlabru)
library(ggplot2)
library(readxl)
library(dplyr)

bru_options_set(control.compute = list(dic = TRUE)) # Activate DIC output

# Read avalanche data and convert date
avalanche_data <- read_excel("~/Desktop/hazard_mapping/Lyngen_detections_attribute.xlsx")
avalanche_data$Dato <- as.Date(avalanche_data$Dato, format = "%Y-%m-%d")

# Define the start and end dates for the desired range
start_date <- as.Date("2020-01-01")
end_date <- as.Date("2020-03-30")

# Subset the dataframe to include data within the specified date range
avalanche_data <- avalanche_data[avalanche_data$Dato >= start_date & avalanche_data$Dato <= end_date, ]
avalanche_data <- avalanche_data[order(avalanche_data$skredTidspunkt), ]

# Read warning data and convert date
warning_data <- read_excel("~/Desktop/hazard_mapping/avalanche-bulletin_v1.xlsx")
warning_data$Date <- as.Date(warning_data$Dato, format = "%Y-%m-%d")

# Prepare point data for exploratory data analysis
point_dta <- data.frame(X = avalanche_data$X, Y = avalanche_data$Y, Date = avalanche_data$Dato)
unique_dates <- unique(point_dta$Date)

# Function to specify decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))

# Load elevation, slope, and aspect data
hoyde_data <- read_excel("~/Desktop/hazard_mapping/hoyde100m.xlsx")
slope_data <- read_excel("~/Desktop/hazard_mapping/slope100m.xlsx")
aspect_data <- read_excel("~/Desktop/hazard_mapping/aspect100m.xlsx")

# Prepare elevation, slope, and aspect data frames
point_dta$time <- as.integer(difftime(point_dta$Date, min(point_dta$Date), units = "days")) + 1

elevData <- data.frame(X = hoyde_data$X, Y = hoyde_data$Y, elev = round(hoyde_data$LYNGEN100M_BAND_ / 1000, 2))
slopeData <- data.frame(X = slope_data$X, Y = slope_data$Y, slope = specify_decimal(slope_data$Slope_lyngen100m_Band_1, 1))
aspectData <- data.frame(X = aspect_data$X, Y = aspect_data$Y, aspect = specify_decimal(aspect_data$Aspect_lyngen100m_Band_1, 1))

# Prepare factor data frames
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

# Apply conversion functions to aspect and slope data
aspectData.Factor$aspect <- sapply(aspectData.Factor$aspect, convert_to_compass_direction_2)
slopeData.Factor$slope <- sapply(slopeData.Factor$slope, convert_grad)
aspectData.Factor$aspect <- factor(aspectData.Factor$aspect)
slopeData.Factor$slope <- factor(slopeData.Factor$slope)

# Convert factor data frames to spatial objects
coordinates(aspectData.Factor) <- ~X + Y
gridded(aspectData.Factor) <- TRUE
proj4string(aspectData.Factor) <- CRS("+proj=utm +zone=33 +datum=WGS84")

coordinates(slopeData.Factor) <- ~X + Y
gridded(slopeData.Factor) <- TRUE
proj4string(slopeData.Factor) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Convert point data to spatial object
coordinates(point_dta) <- c("X", "Y")
proj4string(point_dta) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Define the corners for the window
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

# Create a Polygon object
p <- Polygon(utm_coords)
ps <- Polygons(list(p), "boundary_id")
boundary <- SpatialPolygons(list(ps))
proj4string(boundary) <- CRS("+proj=utm +zone=33 +datum=WGS84")

# Prepare spatiotemporal covariate data
observed_avalanchedata <- avalanche_data
observed_avalanchedata$time <- as.integer(difftime(avalanche_data$Dato, min("2020-01-01"), units = "days")) + 1
warning_data$time <- as.integer(difftime(warning_data$Dato, min("2020-01-01"), units = "days")) + 1

covar <- data.frame(X = aspect_data$X, Y = aspect_data$Y)

# Repeat each row for time points (30 days for January and 90 days for 3 months)
time_points <- rep(1:30, each = nrow(covar))
covar_spacetime <- covar[rep(seq_len(nrow(covar)), times = 30), ]
covar_spacetime$time <- time_points

# Merge covariate data with warning data based on time
warning_data_selected <- warning_data %>% select(time, Utbredelse_SP1, Skredstorrelse_SP1, Utlosbarhet_SP1, Faregrad)
merged_data <- covar_spacetime %>% left_join(warning_data_selected, by = "time")

covar_dta <- merged_data
coordinates(covar_dta) <- ~X + Y
gridded(covar_dta) <- TRUE
proj4string(covar_dta) <- CRS("+proj=utm +zone=33 +datum=WGS84")
covar_dta$Utbredelse_SP1 <- factor(covar_dta$Utbredelse_SP1)
covar_dta$Skredstorrelse_SP1 <- factor(covar_dta$Skredstorrelse_SP1)
covar_dta$Utlosbarhet_SP1 <- factor(covar_dta$Utlosbarhet_SP1)
covar_dta$Faregrad <- factor(covar_dta$Faregrad)

# Functions to extract explanatory variables
f.elev <- function(where) {
  v <- eval_spatial(elevData, where, layer = "elev")
  if (any(is.na(v))) v <- bru_fill_missing(elevData, where, v)
  return(v)
}

f.slope.Factor <- function(where) {
  v <- eval_spatial(slopeData.Factor, where, layer = "slope")
  if (any(is.na(v))) v <- bru_fill_missing(slopeData.Factor, where, v)
  return(v)
}

f.aspect.Factor <- function(where) {
  v <- eval_spatial(aspectData.Factor, where, layer = "aspect")
  if (any(is.na(v))) v <- bru_fill_missing(aspectData.Factor, where, v)
  return(v)
}

f.slope <- function(where) {
  v <- eval_spatial(slopeData, where, layer = "slope")
  if (any(is.na(v))) v <- bru_fill_missing(slopeData, where, v)
  return(v)
}

f.aspect <- function(where) {
  v <- eval_spatial(aspectData, where, layer = "aspect")
  if (any(is.na(v))) v <- bru_fill_missing(aspectData, where, v)
  return(v)
}

# Functions to extract covariate data for stability, frequency, size, and danger
f.cov.stab <- function(X, Y, time) {
  spp <- data.frame(X = X, Y = Y)
  coordinates(spp) <- ~X + Y
  proj4string(spp) <- CRS("+proj=utm +zone=33 +datum=WGS84")
  v <- over(spp, covar_dta[which(covar_dta$time == time), ])
  return(v$Utlosbarhet_SP1)
}

f.cov.freq <- function(X, Y, time) {
  spp <- data.frame(X = X, Y = Y)
  coordinates(spp) <- ~X + Y
  proj4string(spp) <- CRS("+proj=utm +zone=33 +datum=WGS84")
  v <- over(spp, covar_dta[which(covar_dta$time == time), ])
  return(v$Utbredelse_SP1)
}

f.cov.size <- function(X, Y, time) {
  spp <- data.frame(X = X, Y = Y)
  coordinates(spp) <- ~X + Y
  proj4string(spp) <- CRS("+proj=utm +zone=33 +datum=WGS84")
  v <- over(spp, covar_dta[which(covar_dta$time == time), ])
  return(v$Skredstorrelse_SP1)
}

f.cov.danger <- function(X, Y, time) {
  spp <- data.frame(X = X, Y = Y)
  coordinates(spp) <- ~X + Y
  proj4string(spp) <- CRS("+proj=utm +zone=33 +datum=WGS84")
  v <- over(spp, covar_dta[which(covar_dta$time == time), ])
  return(v$Faregrad)
}

# Store the projection information
boundary_proj <- proj4string(boundary)

# Create the mesh for spatial modeling
max.edge <- diff(range(coordinates(point_dta)[, 1])) / (3 * 5)
bound.outer <- diff(range(coordinates(point_dta)[, 1])) / 3

mesh <- inla.mesh.2d(
  loc = coordinates(point_dta),
  boundary = boundary,
  max.edge = c(1, 5) * max.edge,
  offset = c(max.edge, bound.outer),
  cutoff = max.edge / 5,
  crs = CRS(boundary_proj)
)

matern <- inla.spde2.pcmatern(
  mesh,
  prior.sigma = c(0.1, 0.01),
  prior.range = c(10, 0.1)
)

# Define and fit the models
cmp4.1D <- coordinates + time ~ Intercept(1) + danger(f.cov.danger(X, Y, time), model = "factor_contrast") +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = FALSE, hyper = list(theta = list(initial = log(25), fixed = TRUE)))

fit4.1D <- lgcp(cmp4.1D,
                data = point_dta,
                samplers = boundary,
                domain = list(
                  coordinates = mesh,
                  time = seq_len(30)
                )
)


cmp4.1 <- coordinates + time ~ Intercept(1) + stability(f.cov.stab(X, Y, time), model = "factor_contrast") +
  frequency(f.cov.freq(X, Y, time), model = "factor_contrast") +
  size(f.cov.size(X, Y, time), model = "factor_contrast") +
  slope(f.slope.Factor(.data.), model = "factor_contrast") +
  elev(f.elev(.data.), model = "rw1", cyclic = FALSE, hyper = list(theta = list(initial = log(25), fixed = TRUE)))

fit4.1 <- lgcp(cmp4.1,
               data = point_dta,
               samplers = boundary,
               domain = list(
                 coordinates = mesh,
                 time = seq_len(30)
               )
)
