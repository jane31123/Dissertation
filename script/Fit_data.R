# Load the package
library(terra); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(rstudioapi); library(tidyr); library(stringr); library(mgcv); library(gratia)
library(caret); library(ranger); library(pROC); library(lubridate); library(ggforce)
library(geobr); library(viridis)

# set working directory 
setwd("C:/Users/User/OneDrive - University College London/Research Project/Data/arbos_2025")

# dengue surveillance data
# format dates 
dd = read.csv("./dd.csv") 
cd = read.csv("./cd.csv")  
zd = read.csv("./zd.csv")
pop = read.csv("./pop_2016-2024 by municipality.csv")

dc = read.csv("./dc.csv") 
cc = read.csv("./cc.csv")  
zc = read.csv("./zc.csv")

dc<-dd |> filter(Region=="Centro Oeste")
cc<-cd |> filter(Region=="Centro Oeste")
zc<-zd |> filter(Region=="Centro Oeste")  

write.csv(pop, "pop_17-24.csv", row.names=FALSE)
write.csv(cc, "cc.csv", row.names=FALSE)
write.csv(zc, "zc.csv", row.names=FALSE)

# Area
area = read.csv("./area_2024.csv")  
all_dengue <- all_dengue %>%
  mutate(IBGE6 = as.character(IBGE6))
d4 <- d4 %>%
  mutate(IBGE6 = as.character(IBGE6))
c4 <- c4 %>%
  mutate(IBGE6 = as.character(IBGE6))
zd <- zd %>%
  mutate(IBGE6 = as.character(IBGE6))

write.csv(area, "area_2024.csv", row.names=FALSE)

# Create new column, area
dd<-dd%>%
  left_join(area, by="IBGE6")
cd<-cd%>%
  left_join(area, by="IBGE6")
zd<-zd%>%
  left_join(area, by="IBGE6")

# calculate population density of each municipality
dd<-dd%>%mutate(pd=Population/Area)
cd<-cd%>%mutate(pd=Population/Area)
zd<-zd%>%mutate(pd=Population/Area)

# an accompanying shapefile of Brazil districts
shp_brazil = sf::st_read("C:/Users/User/OneDrive - University College London/Research Project/Data/arbos_2025/Brazil_shp_harm_2022.shp")
ggplot() +
  geom_sf(data=shp_brazil, color=NA, fill="grey80")
plot(shp_brazil)

# visualise a histogram of incidence
# what do you notice about this distribution?
summary(dd$Incidence)
summary(cd$Incidence)
summary(zd$Incidence)

# Calculate State level population and incidence rate.
# The one shown originally were municipality-level population and incidence rate.
dd$Date <- as.Date(dd$Date)
cd$Date <- as.Date(cd$Date)
zd$Date <- as.Date(zd$Date)

dd %>%
  dplyr::group_by(Region, State, Date) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE), 
                   Population = sum(Population, na.rm=TRUE),
                   Incidence = NumCases / (Population/100000)) %>%
  ungroup()%>%
  ggplot() + 
  geom_line(aes(Date, Incidence, group=State, color=factor(State))) + 
  scale_x_date(date_labels = "%Y") +
  theme_minimal() + 
  facet_wrap_paginate(~Region, nrow = 1, ncol = 1, page = 5) + # plots separate subplot for each state 
  xlab("Year") + 
  ylab("Dengue incidence per 100,000")+
  labs(color = "State") +  # Change legend title
  ggtitle("Dengue State level Incidence Rate divided by 5 Regions")  # Add plot title

ggsave(filename="d53.jpg", device="jpg", units="in", width=10, height=5, dpi=600)

cd %>%
  dplyr::group_by(Region, State, Date) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE), 
                   Population = sum(Population, na.rm=TRUE),
                   Incidence = NumCases / (Population/100000)) %>%
  ungroup()%>%
  ggplot() + 
  geom_line(aes(Date, Incidence, group=State, color=factor(State))) + 
  scale_x_date(date_labels = "%Y") +
  theme_minimal() + 
  facet_wrap_paginate(~Region, nrow = 1, ncol = 2, page = 2) + # plots separate subplot for each state 
  theme(legend.position = "none")+
  xlab("Year") + 
  ylab("Chikungunya incidence per 100,000")+ 
  ggtitle("Chikungunya State level Incidence Rate divided by 5 Regions")  # Add plot title

ggsave(filename="C52.jpg", device="jpg", units="in", width=10, height=5, dpi=600)

zd %>%
  dplyr::group_by(Region, State, Date) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE), 
                   Population = sum(Population, na.rm=TRUE),
                   Incidence = NumCases / (Population/100000)) %>%
  ungroup()%>%
  ggplot() + 
  geom_line(aes(Date, Incidence, group=State, color=factor(State))) + 
  scale_x_date(date_labels = "%Y") +
  theme_minimal() + 
  facet_wrap_paginate(~Region, nrow = 1, ncol = 1, page = 5) + # plots separate subplot for each state 
  xlab("Year") + 
  ylab("Zika incidence per 100,000")+
  labs(color = "State") +  
  ggtitle("Zika State level Incidence Rate divided by 5 Regions")  # Add plot title
ggsave(filename="Z53.jpg", device="jpg", units="in", width=10, height=5, dpi=600)

#-------------------------------------------------------------------------------
# visualize 2017-2024 Aedes-borne disease cumulative incidence rate by municipality
# level
table <- table %>%
  mutate(IBGE6 = as.character(IBGE6))

shp_brazil<- shp_brazil%>%
  mutate(IBGE6 = as.character(IBGE6))

# Plot incidence level on map
# add up total cases of three disease in 8 years by municipality
total = all %>%
  dplyr::group_by(IBGE6) %>%
  dplyr::summarise(total_case = sum(NumCases))%>%
  ungroup() 

# extract 2017 population of each municipality 
pop = read.csv("./pop_2016-2024 by municipality.csv") 
pop_17<-pop |> filter(Year=="2017")

# make the dataframe contain shape file
table<-total |> left_join(pop_17, by= "IBGE6")   
table<-table |> mutate(incidence=(total_case/Population)*100000) 
shp_brazil<-sf::st_transform(shp_brazil, crs = 4326)

# combine with shapefile
shp_total = shp_brazil %>%
  dplyr::left_join(table, by="IBGE6") 

# exclude 2 rows of NAs
shp_total_cleaned <- shp_total |>  drop_na()

# plot with Brazil country as background
# Original plot with the region boundaries added on top
ggplot() +
  geom_sf(data = shp_total_cleaned, aes(fill = incidence), color = "grey", size = 0.1) +
  geom_sf(data = region, fill = NA, color = "black", size = 1.2) + # This adds the region borders
  geom_sf_text(data = region, aes(label = name_region), color = "black", size = 3.5) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Aedes-borne disease\nincidence") +
  theme_minimal() +
  ggtitle("Cumulative incidence rate, 2017-2024") 

ggsave(filename="incidence rate of 5 regions.jpg", device="jpg", units="in", width=7, height=5, dpi=600)

library(ggplot2)
library(sf)
library(viridis)

dissertation_incidence_map <- ggplot() +
  # Municipal polygons with incidence data
  geom_sf(data = shp_total_cleaned, aes(fill = incidence), 
          colour = "white", size = 0.05) +  # Thinner, white boundaries
  
  # Region boundaries (more prominent)
  geom_sf(data = region, fill = NA, 
          colour = "black", size = 0.8, linetype = "solid") +
  
  # Region labels with better positioning
  geom_sf_text(data = region, aes(label = name_region), 
               colour = "black", size = 3.2, fontface = "bold",
               check_overlap = TRUE) +
  
  # Professional colour scale
  scale_fill_viridis_c(
    option = "inferno", 
    direction = -1,
    name = "Incidence rate\n(per 100,000)",  # More specific
    trans = "sqrt",  # Square root transformation for better visual distribution
    labels = scales::comma_format(accuracy = 0.1),
    guide = guide_colorbar(
      barwidth = 12,
      barheight = 0.8,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  # Clean coordinate system
  coord_sf(expand = FALSE, crs = st_crs(4326)) +
  
  # Academic-style labels
  labs(
    title = "Cumulative Aedes-borne Disease Incidence, 2017–2024",
    subtitle = "Incidence rate per 100,000 population by municipality",
    caption = "Source: Brazilian Ministry of Health (DATASUS); IBGE administrative boundaries"
  ) +
  
  # Clean academic theme
  theme_void() +
  theme(
    # Title formatting
    plot.title = element_text(size = 14, face = "bold", 
                              hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, 
                                 margin = margin(b = 15), colour = "grey40"),
    
    # Legend formatting
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    legend.box.margin = margin(t = 15),
    legend.margin = margin(t = 5),
    
    # Caption
    plot.caption = element_text(size = 8, hjust = 0, 
                                margin = margin(t = 10), colour = "grey50"),
    
    # Margins and background
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Remove any remaining grid
    panel.grid = element_blank()
  )

print(dissertation_incidence_map)

# Replace scale_fill_viridis_c() with:
scale_fill_distiller(
  palette = "Spectral",
  direction = -1,
  name = "Incidence rate\n(per 100,000)",
  trans = "sqrt",
  labels = scales::comma_format(accuracy = 0.1)
)

# For better accessibility and printing
scale_fill_viridis_c(
  option = "plasma",  # Better for colorblind readers
  direction = -1,
  name = "Incidence rate\n(per 100,000)",
  trans = "sqrt"
)

# High-resolution export
ggsave("incidence_map_dissertation.png", dissertation_incidence_map,
       width = 10, height = 8, dpi = 300, bg = "white")

# Vector format for LaTeX
ggsave("incidence_map_dissertation.pdf", dissertation_incidence_map,
       width = 10, height = 8, device = "pdf")

# For presentations
ggsave("incidence_map_presentation.png", dissertation_incidence_map,
       width = 12, height = 9, dpi = 150, bg = "white")


