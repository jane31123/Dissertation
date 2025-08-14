# Load the package
library(terra); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(rstudioapi); library(tidyr); library(stringr); library(mgcv); library(gratia)
library(caret); library(ranger); library(pROC); library(lubridate); library(ggforce)

# set working directory 
setwd("C:/Users/User/OneDrive - University College London/Research Project/Data/arbos_2025")

# load zika dataframe 
zc= read.csv("./zc.csv") 

zc$Month <- as.factor(zc$Month)
zc$IBGE6 <- as.factor(zc$IBGE6)
zc$Year  <- as.factor(zc$Year) 
zc$Date <-as.Date(zc$Date)

# subset to only high-incidence locations
# > 80th percentile of incidence
hi = zc %>%
  dplyr::group_by(IBGE6) %>%
  dplyr::summarise(Incidence = mean(Incidence)) %>%
  dplyr::arrange(desc(Incidence))

# identify 80th percentile
pc80 = quantile(hi$Incidence, 0.8)

# filter down to just the high incidence municipalities
hi = hi %>% dplyr::filter(Incidence > pc80)
zzc = zc %>%
  dplyr::filter(IBGE6 %in% hi$IBGE6)

# format date and set IBGE6 and year as factor
zzc$Date  =  as.Date(zzc$Date)
zzc$IBGE6 <- as.factor(zzc$IBGE6)
zzc$Year  <- as.factor(zzc$Year) 
zzc$State =  as.factor(zzc$State)
zzc$Month =  as.numeric(zzc$Month)

# ------------------plot zika incidence - just from 3 states------------
nzc<-zc |> dplyr::select(IBGE6, Municipality, State, Population, NumCases, Year, Date, Month)

# Calculate state-level monthly aggregates and incidence rate
state_monthly <- nzc %>%
  group_by(State, Date) %>%
  summarise(
    total_population = sum(Population, na.rm = TRUE),
    total_cases = sum(NumCases, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    incidence = (total_cases / total_population) * 100000
  )

# Verify the aggregation worked correctly
state_monthly %>%
  filter(State == "Mato Grosso do Sul", Date == "2017-01-15") 

# Compare with original municipal data for the same period
nzc %>%
  filter(State == "Mato Grosso do Sul", Year == 2017, Month == 1) %>%
  summarise(
    check_population = sum(Population),
    check_cases = sum(NumCases)
  )

pz <- ggplot(state_monthly, aes(x = Date, y = incidence, colour = State)) +
  geom_line(size = 1.2, alpha = 0.8) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = c(0.02, 0)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(title = "Zika Incidence Rate by State (2017-2024)",
       x = "Year",
       y = "Incidence Rate \n(per 100,000 population)",
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

print(pz) 

 
#==============================================================================
# create offset for log population (100,000s) so you are modelling incidence
zzc$offset = log(zzc$Population/100000)

# check number of zeroes - not zero inflated now
table(zzc$NumCases == 0)


################################################

# 1. fit "baseline" model with no climate info

# includes offset, monthly thin plate spline stratified by State, year*state interaction, random effect of IBGE6
# this function also checks fitting time
s = Sys.time()
m_base <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re"),
              family = "nb",
              method = "REML",
              data = zzc)
e = Sys.time()


# on my computer this takes 10 seconds - much better!
e-s 

#-------------------------------------------------------------------------------
# test lags for temperature
# adding 1 climate variable at a time to the baseline model

m_t0 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t5 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_t6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)



# -------------------------------------precipitation----------------------------------
m_p0 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p5 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_p6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

#-------------------------------SPI-6----------------------------------
m_s0 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp")+ s(spi6, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s5 <- gam(NumCases ~ -1 + offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

m_s6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_6m, bs="tp") + s(spi6_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = zzc)

# compare AIC within temperature models, tmean gets best result 
AIC(m_base, m_t0, m_t1, m_t2, m_t3, m_t4, m_t5, m_t6, 
    m_p0, m_p1, m_p2, m_p3, m_p4, m_p5, m_p6, m_s0, m_s1, m_s2, m_s3, m_s4, m_s5, m_s6)

# extract and plot spline from best fitting model
# n.b. apply exponential transform to gain incidence rate ratio
draw(m_s5, fun=exp, select="s(tmean_6m)", caption=FALSE)+ 
  labs(title="(A) Association between lagged Temperature and Disease Incidence Rate", 
       x="Monthly Mean Temperature with 6-month lag (°C)", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="zika-T6.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s5, fun=exp, select="s(pmonth)", caption=FALSE)+ 
  labs(title="(B) Association between Precipitation and Disease Incidence Rate", 
       x="Monthly Precipitation (mm)", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="zika-P0.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s5, fun=exp, select="s(spi6_5m)", caption=FALSE)+ 
  labs(title="(C) Association between SPI-6 and Disease Incidence Rate", 
       x="SPI-6 with 5 months lag effect", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="zika-S5.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s5, fun=exp, select=c("s(Month):StateMato Grosso") , caption=FALSE) + 
  labs(title="Seasonal Variation of Zika Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="MG.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s5, fun=exp, select=c("s(Month):StateGoiás") , caption=FALSE) + 
  labs(title="Seasonal Variation of Zika Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="GO.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s5, fun=exp, select=c("s(Month):StateMato Grosso do Sul") , caption=FALSE) + 
  labs(title="Seasonal Variation of Zika Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="MS.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

#========================================================================================
# Example: assign rough population sizes (or use actual means from your data)
nn<-zc |> dplyr::select(IBGE6, State, Region, Municipality, Population, Year, Month)



# Extract model predictions for year-month combinations
new_data <- expand.grid(
  Month = 1:12,
  Year = 2017:2024,
  State = c("Goiás", "Mato Grosso", "Mato Grosso do Sul"),
  IBGE6 = "dummy", # placeholder
  tmean_6m = mean(zzc$tmean_6m),
  pmonth = mean(zzc$pmonth),
  spi6_5m = mean(zzc$spi6_5m)
)

# create state-year population table
state_pop <- nn %>%
  filter(Month==1) %>%
  group_by(State, Year) %>%
  summarise(State_Pop = sum(Population)) |> 
  ungroup()

print(state_pop)

# Merge with new_data by State and Year
new_data <- merge(new_data, state_pop, by = c("State", "Year"))

# Create the Offset
new_data$offset <- log(new_data$State_Pop / 100000)

# Predict and plot
predictions <- predict(m_s5, newdata=new_data, type="response")
new_data$predicted <- predictions

ggplot(new_data, aes(x=Month, y=Year, fill=predicted)) +
  geom_tile() +
  facet_wrap(~State) +
  scale_fill_viridis_c(name="Fitted\nIncidence rate") +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=seq(2017,2024,1)) +
  theme_minimal() +
  labs(title="Fitted zika incidence by month and year across states")

# ------------------------------Diagnosis-------------------------------------
# check the output
summary(m_s5) 

# check residuals
par(mfrow = c(2, 2))  # Set up 2 rows and 2 columns
gam.check(m_s5)

 
save(list = ls(all.names = TRUE), file = "zika.RData")  
