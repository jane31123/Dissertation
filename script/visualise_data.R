# or use the development version with latest features
utils::remove.packages('geobr')
remotes::install_github("ipeaGIT/geobr", subdir = "r-package")

# Packages dependency
library(terra); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(rstudioapi); library(tidyr); library(stringr); library(mgcv); library(gratia)
library(raster);library(lubridate); library(purrr); library(rhdf5); library(exactextractr)
library(geobr); library(scales)

# ----------- data read in and visualisation ------------

# set working directory
setwd("C:/Users/User/OneDrive - University College London/Research Project/Data/arbos_2025")

bra =read.csv("C:/Users/User/OneDrive - University College London/Research Project/Data/brazil_climvars_era5.csv")
ds = read.csv("./ds.csv")
cs = read.csv("./cs.csv")
zs = read.csv("./zs.csv")
climate = read.csv("./climate_2015-2024.csv")
climate<- climate |> dplyr::select(-spi1,-spi1_1m, -spi1_2m,-spi1_3m)
write.csv(climate, "climate_2015-2024.csv", row.names=FALSE)

# ----------- for dengue -------------

# dengue data - monthly case counts for all municipalities
# IBGE6 = unique municipality code
dd = read.csv("./dengue_monthlyTS_2001_2024.csv") %>%
  dplyr::mutate(Date = as.Date(Date)) 

# use geographic lookup file to get name, state and region of municipality
geo = read.csv("./geo_lookup_table.csv")
dd = dplyr::left_join(dd, geo) 

# add populations and calculate incidence 
pop = read.csv("./pop_2016-2024 by municipality.csv")

dd <- dd %>%   # Remove existing Population column
  dplyr::left_join(pop, by = c("IBGE6", "Year"))%>%
  dplyr::filter(Year%in% 2016:2024)

# calculate incidence per 100,000 persons
dd$Incidence = dd$NumCases / (dd$Population/100000)

# ----------- same as above, for Chikungunya -------------

# Chikungungya data - monthly case counts for all municipalities
# IBGE6 = unique municipality code
cd = read.csv("./chikungunya_monthlyTS_2017_2024.csv") %>%
  dplyr::mutate(Date = as.Date(Date))

# use geographic lookup file to get name, state and region of municipality
cd = dplyr::left_join(cd, geo)

# add populations and calculate incidence 
cd = dplyr::left_join(cd, pop)

# calculate incidence per 100,000 persons
cd$Incidence = cd$NumCases / (cd$Population/100000) 

# summarise totals over entire monitoring period and map
totals_cd = cd %>%
  dplyr::group_by(IBGE6) %>%
  dplyr::summarise(Total = sum(NumCases, na.rm=TRUE))%>%
  ungroup()


# ----------- same as above, for Zika -------------

# Zika data - monthly case counts for all municipalities
# IBGE6 = unique municipality code
zd = read.csv("./zika_monthlyTS_2016_2024.csv") %>%
  dplyr::mutate(Date = as.Date(Date))

# use geographic lookup file to get name, state and region of municipality
geo = read.csv("./geo_lookup_table.csv")
zd = dplyr::left_join(zd, geo)

# add populations and calculate incidence 
pop = read.csv("./pop.csv")
zd = dplyr::left_join(zd, pop)

# calculate incidence per 100,000 persons
zd$Incidence = zd$NumCases / (zd$Population/100000)

# ------------ shapefile ------------

# brazil shapefile; note the coordinate ref system
shp = sf::st_read("./Brazil_shp_harm_2022.shp")
  sf::st_transform(shp, crs = 4326) 

#---------------------------------------geobr-----------------------------------
# Available data sets
datasets <- list_geobr()
datasets

# View Municipality
muni <- read_municipality(year = 2022,
                          showProgress = FALSE)
head(muni)
nrow(muni)
plot(muni["code_muni"])

muni <- muni %>%
  mutate(IBGE6 = substr(as.character(code_muni), 1, 6))

# Read 27 States
state <- read_state(code_state="all", year=2020)

# Read 5 Regions
region<- read_region(year=2020) 

# Filter the data to get only the "Centro Oeste" region
# We assume your 'state' sf object is loaded in your environment.
co <- state %>%
  filter(name_region == "Centro Oeste")

muni |> sf::st_transform(muni, crs = 4326)
region |> sf::st_transform(region, crs = 4326)
state |> sf::st_transform(state, crs = 4326)

# Define the custom color palette
# Create a named vector where the names match the 'name_state' column. 
color_palette <- c(
  "Distrito Federal" = "coral",
  "Goiás" = "lightgreen",
  "Mato Grosso" = "lightblue",
  "Mato Grosso Do Sul" = "mediumpurple"
)


# Create publication-ready map
dissertation_map <- ggplot() +
  # State polygons with clean boundaries
  geom_sf(data = co, aes(fill = name_state), 
          colour = "white", linewidth = 0.3) +
  
  # Municipal boundaries (subtle)
  geom_sf(data = co_muni, fill = NA, 
          colour = "grey70", linewidth = 0.1, alpha = 0.7) +
  
  # Professional colour palette
  scale_fill_manual(values = color_palette) +
  
  # Coordinate system with proper projection
  coord_sf(crs = st_crs(4326), 
           expand = FALSE,  # Remove extra white space
           datum = sf::st_crs(4326)) +
  
  # Professional labels
  labs(
    title = "Centro-Oeste Region of Brazil",
    subtitle = "States and Federal District with municipal boundaries",
    fill = "State/Federal District",
    caption = "Source: Brazilian Institute of Geography and Statistics (IBGE)"
  ) +
  
  # Clean academic theme
  theme_minimal() +
  theme(
    # Title formatting
    plot.title = element_text(size = 14, face = "bold", 
                              hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, 
                                 margin = margin(b = 15), colour = "grey40"),
    
    # Legend formatting
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box.margin = margin(t = 15),
    legend.key.size = unit(0.8, "cm"),
    
    # Caption formatting
    plot.caption = element_text(size = 9, hjust = 0, 
                                margin = margin(t = 10), colour = "grey50"),
    
    # Overall plot margins
    plot.margin = margin(20, 20, 20, 20),
    
    # Background
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(dissertation_map)
ggsave("co.png", dissertation_map,
       width = 8, height = 6, dpi = 300, bg = "white")
 


# ----------- climate --------------
# monthly temperature data
# read in and rename each raster layer to format "X-DATE" 
# convert from kelvin to celcius
temp <- rast("./temperature.grib")-273.15 
 
plot(temp[[1:6]]) 
terra::crs(temp) <- "EPSG:4326"
terra::crs(prec) <- "EPSG:4326" 
state<-sf::st_transform(state, crs = 4326)
co<-state |> filter(name_region=="Centro Oeste")

# Crops the raster temp to the spatial extent of the reprojected shapefile xx
temp_crop <- terra::crop(temp, co) 
temp_crop[[1]]
names(temp) = paste("X", substr(terra::time(temp), 1, 10), sep="")
names(temp_crop) = paste("X", substr(terra::time(temp_crop), 1, 10), sep="")

# TEMPERATURE
# extract mean climate values for each polygon and location
ext_temp <- terra::extract(temp_crop, co, fun = mean, na.rm = TRUE, ID=FALSE) 
ext_temp$code_state = co$code_state # add unique ID column 

# Move IBGE6 to first row
ext_temp <- ext_temp %>%
  dplyr::relocate(code_state)

# transform to long format
ext_temp_long <- ext_temp %>%
  tidyr::pivot_longer(
    cols = -code_state,
    names_to = "date",
    values_to = "T"
  ) %>%
  mutate(
    # Tell as.Date the specific format of your date strings
    date = as.Date(gsub("^X", "", date), format = "%Y-%m-%d") 
  )

# add month and year column
ext_temp_long<-ext_temp_long |> dplyr::mutate(year = lubridate::year(date), 
                                                  month = lubridate::month(date))

# Add abbrev_state and name_state columns
ext_temp_long <- ext_temp_long %>%
  left_join(co %>% dplyr::select(code_state, abbrev_state, name_state), 
            by = "code_state")


ext_temp_long<-ext_temp_long |> filter(year%in%c(2017:2024))

# plot state level temperature seasonal change through 8 years
p <- ggplot(ext_temp_long, aes(x = date, y = T, colour = name_state)) +
  geom_line(size = 1.2, alpha = 0.8) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = c(0.02, 0)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(title = "Monthly mean Temperature by State (2017-2024)",
       x = "Year",
       y = "Monthly mean temperature (°C)",
       colour = "State") +  # UK English spelling
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    panel.grid.minor = element_blank() 
  ) +
  guides(colour = guide_legend(nrow = 1))

print(p)



# --------------------------Precipitation---------------------------------------
# monthly precip data and rename layers
# load GRIB data
prec <- rast("./precipitation.grib")*1000

# label column names as their recorded time 
names(prec) = paste("X", substr(terra::time(prec), 1, 10), sep="") 
names(prec_crop) = paste("X", substr(terra::time(prec_crop), 1, 10), sep="")
prec_crop[[1]] 

# transform CRS
terra::crs(prec) <- "EPSG:4326"

# Crops the raster temp to the spatial extent of the reprojected shapefile xx
prec_crop <- terra::crop(prec, co)
prec_crop[[1]] 

# extract mean climate values for each polygon and location
ext_prec <- terra::extract(prec_crop, co, fun = mean, na.rm = TRUE, ID=FALSE) 
ext_prec$code_state = co$code_state # add unique ID column 

# Move IBGE6 to first row
ext_prec <- ext_prec %>%                    
  dplyr::select(code_state, everything())

# transform to long format
ext_prec_long <- ext_prec %>%
  tidyr::pivot_longer(
    cols = -code_state,                         # all columns *except* IBGE6
    names_to = "date",
    values_to = "P"
  ) %>%
  mutate(
    date = as.Date(gsub("^X", "", date))   # remove the "X" prefix from column names
  )

# transform daily precipitation to monthly precipitation
ext_prec_long <- ext_prec_long %>%
  mutate(
    days_in_month = days_in_month(date),
    pmonth = P * days_in_month
  )

# create month and year column
ext_prec_long<-ext_prec_long |> dplyr::mutate(year = lubridate::year(date), 
              month = lubridate::month(date))
 
ext_prec_long<-ext_prec_long |> filter(year%in%c(2017:2024))

# Add abbrev_state and name_state columns
ext_prec_long <- ext_prec_long %>%
  left_join(co %>% dplyr::select(code_state, abbrev_state, name_state), 
            by = "code_state")

# plot state level temperature seasonal change through 8 years
pp <- ggplot(ext_prec_long, aes(x = date, y = pmonth, colour = name_state)) +
  geom_line(size = 1.2, alpha = 0.8) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = c(0.02, 0)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(title = "Monthly Precipitation by State (2017-2024)",
       x = "Year",
       y = "Monthly Precipitation (mm)",
       colour = "State") +  # UK English spelling
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    panel.grid.minor = element_blank() 
  ) +
  guides(colour = guide_legend(nrow = 1))

print(pp)

# create climate table
climate = dplyr::left_join(ext_temp_long, ext_precip_long) %>%
  dplyr::mutate(year = lubridate::year(date), 
                month = lubridate::month(date))
climate <- climate %>%
  mutate(IBGE6 = as.character(IBGE6))
write.csv(climate, "climate_t+p.csv", row.names=FALSE)

#-------------------------------------------------------------------------------
cd <- cd %>%
  mutate(IBGE6 = as.character(IBGE6)) 

dd <- dd %>%
  mutate(IBGE6 = as.character(IBGE6))

zd <- zd %>%
  mutate(IBGE6 = as.character(IBGE6))

temp_lag3 <- temp_lag3 %>%
  mutate(IBGE6 = as.character(IBGE6))

#-------------------------------------------------------------------------------
# set working directory
setwd("C:/Users/User/OneDrive - University College London/Research Project/Data/SPI") 

# Load the NetCDF file
nc1 <- terra::rast("spa01_m_gdo_20010101_20011221_t.nc")
nc24<- terra::rast("spa01_m_gdo_20240101_20241221_t.nc")

# View metadata
nc1
nc24

# List all raster files (adjust the path and pattern if needed)
nc_files <- list.files(path = "C:/Users/User/OneDrive - University College London/Research Project/Data/SPI", 
                       pattern = "\\.nc$", full.names = TRUE)

# Define extent: xmin, xmax, ymin, ymax
brazil_extent <- ext(-73.95, -28.85, -33.75, 5.25)

# Crop each NetCDF file and store results
cropped_list <- lapply(nc_files, function(file) {
  r <- rast(file)
  crop(r, brazil_extent)
})

# Assign file names as names of list elements
names(cropped_list) <- basename(nc_files)

# Split into SPI-1 and SPI-6 groups
spi1_list <- cropped_list[grepl("spa01", names(cropped_list))]
spi6_list <- cropped_list[grepl("spa06", names(cropped_list))]

names(spi1_list)
names(spi6_list)

# Rename Layers Based on Time Values
spi1_list <- lapply(spi1_list, function(r) {
  time_vec <- terra::time(r)
  new_names <- paste0("X", format(time_vec, "%Y-%m-%d"))
  names(r) <- new_names
  return(r)
})

spi6_list <- lapply(spi6_list, function(r) {
  time_vec <- terra::time(r)
  new_names <- paste0("X", format(time_vec, "%Y-%m-%d"))
  names(r) <- new_names
  return(r)
})

# Ensure the CRS matches the raster (should be EPSG:4326)
shp <- st_transform(shp, crs = crs(spi1_list[[1]]))

# Function to crop and mask a raster to shp_filtered
cropped_region <- function(raster) {
  terra::mask(terra::crop(raster, shp), shp)
}

# Apply to SPI-1 list
spi1_cropped <- lapply(spi1_list, cropped_region)

# Apply to SPI-6 list
spi6_cropped <- lapply(spi6_list, cropped_region)

# Extract mean values for each raster in the list
spi1_extract <- lapply(spi1_cropped, function(r) {
  terra::extract(r, shp, fun = mean, na.rm = TRUE, ID = FALSE)
})

spi6_extract <- lapply(spi6_cropped, function(r) {
  terra::extract(r, shp, fun = mean, na.rm = TRUE, ID = FALSE)
})
  
# spi_extract is a list of 24 data frames, each with 2261 rows and 12 columns
length(spi1_extract)  # should be 25
nrow(spi1_extract[[1]])  # should be 5570
ncol(spi1_extract[[1]])    
ncol(spi6_extract[[1]])

# Bind columns from all files side-by-side (each file contributes 12 columns)
spi6_wide <- do.call(cbind, spi6_extract)
spi1_wide <- do.call(cbind, spi1_extract)

# Rename all columns
names(spi6_wide) <- str_extract(names(spi6_wide), "X\\d{4}-\\d{2}-\\d{2}")
names(spi1_wide) <- str_extract(names(spi1_wide), "X\\d{4}-\\d{2}-\\d{2}")

# Add IBGE6 column
spi6_wide$IBGE6 = shp$IBGE6
spi6_wide <- spi6_wide[, c("IBGE6", setdiff(names(spi6_wide), "IBGE6"))]

spi1_wide$IBGE6 = shp$IBGE6
spi1_wide <- spi1_wide[, c("IBGE6", setdiff(names(spi1_wide), "IBGE6"))]

spi6_long <- spi6_wide %>%
  tidyr::pivot_longer(
    cols = -IBGE6,             # all columns *except* IBGE6
    names_to = "date",
    values_to = "spi6"
  ) %>%
  mutate(
    date = as.Date(gsub("^X", "", date))   # remove the "X" prefix from column names
  ) 

spi1_long <- spi1_wide %>%
  tidyr::pivot_longer(
    cols = -IBGE6,             # all columns *except* IBGE6
    names_to = "date",
    values_to = "spi1"
  ) %>%
  mutate(
    date = as.Date(gsub("^X", "", date))   # remove the "X" prefix from column names
  )

spi1_long_with_mean <- spi1_long %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  group_by(IBGE6, year_month) %>%
  mutate(spi1_mean = mean(spi1, na.rm = TRUE)) %>%
  ungroup()

spi1_monthly_mean <- spi1_long_with_mean %>%
  dplyr::select(IBGE6, year_month, spi1_mean) %>%
  distinct()

spi1_monthly_mean<-spi1_monthly_mean %>%                      # Remove unwanted columns
  tidyr::separate(year_month, into = c("year", "month"), sep = "-")

spi6_long<-spi6_long %>%                      # Remove unwanted columns
  tidyr::separate(date, into = c("year", "month", "day"), sep = "-")%>%
  dplyr::select(-day)

# Convert year and month columns to double
spi6_long$year <- as.numeric(spi6_long$year)
spi6_long$month <- as.numeric(spi6_long$month)
spi1_monthly_mean$year <- as.numeric(spi1_monthly_mean$year)
spi1_monthly_mean$month <- as.numeric(spi1_monthly_mean$month) 

climate<- climate%>% dplyr::left_join(spi1_monthly_mean, by = c("year", "month", "IBGE6")) %>%
  dplyr::left_join(spi6_long, by = c("year", "month", "IBGE6"))%>%
  rename(spi1=spi1_mean)

climate_subset <- climate %>%
  dplyr::select(IBGE6, year, month, tmean_3m, tmean_4m, tmean_5m, tmean_6m)
zd <- zd |> 
  left_join(climate_subset, by = c("IBGE6", "Year"="year", "Month"="month"))
write.csv(zd, "zd.csv", row.names=FALSE) 

# create a 1 month lag variable
lag1 = climate %>% 
  dplyr::select(IBGE6, date, tmean, pmonth, spi1, spi6) %>%
  dplyr::mutate(date = date %m+% months(1)) %>%
  dplyr::rename(tmean_1m = tmean,
                pmonth_1m = pmonth, 
                spi1_1m = spi1, 
                spi6_1m = spi6)
climate = dplyr::left_join(climate, lag1)

lag2 = climate %>% 
  dplyr::select(IBGE6, date, tmean, pmonth, spi1, spi6) %>%
  dplyr::mutate(date = date %m+% months(2)) %>%
  dplyr::rename(tmean_2m = tmean,
                pmonth_2m = pmonth, 
                spi1_2m = spi1, 
                spi6_2m = spi6)
climate = dplyr::left_join(climate, lag2)

lag3 = climate %>% 
  dplyr::select(IBGE6, date, pmonth, spi1, spi6) %>%
  dplyr::mutate(date = date %m+% months(3)) %>%
  dplyr::rename(pmonth_3m = pmonth, 
                spi1_3m = spi1, 
                spi6_3m = spi6)

temp_lag3 <- climate %>% 
  mutate(date = as.Date(date)) %>%                  # ensure it's Date
  dplyr::select(IBGE6, date, tmean) %>%
  mutate(date = date %m+% months(3)) %>%            # add 3-month lag
  rename(tmean_3m = tmean)

temp_lag4 <- climate %>% 
  dplyr::select(IBGE6, date, tmean) %>%
  mutate(date = date %m+% months(4)) %>%            # add 4-month lag
  rename(tmean_4m = tmean)

temp_lag5 <- climate %>% 
  dplyr::select(IBGE6, date, tmean) %>%
  mutate(date = date %m+% months(5)) %>%            # add 5-month lag
  rename(tmean_5m = tmean)

temp_lag6 <- climate %>% 
  dplyr::select(IBGE6, date, tmean) %>%
  mutate(date = date %m+% months(6)) %>%            # add 6-month lag
  rename(tmean_6m = tmean)

climate = dplyr::left_join(climate, temp_lag6)

p_lag6 = p_lag6 %>%
  dplyr::mutate(date = as.Date(date))
                
p_lag4 <- climate %>% 
  dplyr::select(IBGE6, date, pmonth) %>%
  mutate(date = date %m+% months(4)) %>%            # add 4-month lag
  rename(pmonth_4m = pmonth)

p_lag5 <- climate %>% 
  dplyr::select(IBGE6, date, pmonth) %>%
  mutate(date = date %m+% months(5)) %>%            # add 5-month lag
  rename(pmonth_5m = pmonth)

p_lag6 <- climate %>% 
  dplyr::select(IBGE6, date, pmonth) %>%
  mutate(date = date %m+% months(6)) %>%            # add 6-month lag
  rename(pmonth_6m = pmonth)
climate = dplyr::left_join(climate, p_lag6, by=c("IBGE6", "date"))

lag4 = climate %>% 
  dplyr::select(IBGE6, date, spi6) %>%
  dplyr::mutate(date = date %m+% months(4)) %>%
  dplyr::rename(spi6_4m = spi6)
climate = dplyr::left_join(climate, lag4)

lag5 = climate %>% 
  dplyr::select(IBGE6, date, spi6) %>%
  dplyr::mutate(date = date %m+% months(5)) %>%
  dplyr::rename(spi6_5m = spi6)
climate = dplyr::left_join(climate, lag5)

lag6 = climate %>% 
  dplyr::select(IBGE6, date, spi6) %>%
  dplyr::mutate(date = date %m+% months(6)) %>%
  dplyr::rename(spi6_6m = spi6)
climate = dplyr::left_join(climate, lag6)

# save
write.csv(zc, "zc.csv", row.names=FALSE)
 

cd <- cd %>%
  left_join(
    climate %>% select(IBGE6, year, month, tmean_3m, tmean_4m, tmean_5m, tmean_6m),
    by = c("IBGE6", "Year"="year", "Month"="month") 
  )

cd <- cd %>%
  left_join(
    climate,
    by = c("IBGE6", "Year" = "year", "Month" = "month") 
  )
zd <- zd %>%
  left_join(
    climate,
    by = c("IBGE6", "Year" = "year", "Month" = "month") 
  ) 

save(list = ls(all.names = TRUE), file = "my_workspace.RData")

