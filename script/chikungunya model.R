
# ------------- SUGGESTED MODIFICATION TO CODE FOR SPEED -------------------

# just fit the model to high-incidence municipalities
# identify municipalities with mean incidence above the 80th percentile
# fit model to just these - keeps most of the rich data, excludes most of the zero-inflation
# and speeds up model fitting
# this can be applied to each disease - means you'll fit models to slightly different
# set of locations for each disease 

# Load the package
library(terra); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(rstudioapi); library(tidyr); library(stringr); library(mgcv); library(gratia)
#library(caret); library(ranger); library(pROC); library(lubridate); library(ggforce)

# set working directory 
setwd("C:/Users/User/OneDrive - University College London/Research Project/Data/arbos_2025")

# load chikungunya dataframe
cc = read.csv("./cc.csv")

# subset to only high-incidence locations
# > 80th percentile of incidence
hi = cc %>%
  dplyr::group_by(IBGE6) %>%
  dplyr::summarise(Incidence = mean(Incidence)) %>%
  dplyr::arrange(desc(Incidence))

# identify 80th percentile
pc80 = quantile(hi$Incidence, 0.8)

# filter down to just the high incidence municipalities
hi = hi %>% dplyr::filter(Incidence > pc80)
ccc = cc %>%
  dplyr::filter(IBGE6 %in% hi$IBGE6)

# format date and set IBGE6 and year as factor
ccc$Date  =  as.Date(ccc$Date)
ccc$IBGE6 <- as.factor(ccc$IBGE6)
ccc$Year  <- as.factor(ccc$Year) 
ccc$State =  as.factor(ccc$State)
ccc$Month =  as.numeric(ccc$Month)


# ------------------plot chikungunya incidence - just from 3 states------------
ncc<-cc |> dplyr::select(IBGE6, Municipality, State, Population, NumCases, Year, Date, Month)

# Calculate state-level monthly aggregates and incidence rate
state_monthly <- ncc %>%
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
ncc %>%
  filter(State == "Mato Grosso do Sul", Year == 2017, Month == 1) %>%
  summarise(
    check_population = sum(Population),
    check_cases = sum(NumCases)
  )

pc <- ggplot(state_monthly, aes(x = Date, y = incidence, colour = State)) +
  geom_line(size = 1.2, alpha = 0.8) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = c(0.02, 0)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(title = "Chikungunya Incidence Rate by State (2017-2024)",
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

print(pc) 

#=============================================================================
# create offset for log population (100,000s) so you are modelling incidence
ccc$offset = log(ccc$Population/100000)

# check number of zeroes - not zero inflated now
table(ccc$NumCases == 0)


################################################

# 1. fit "baseline" model with no climate info

# includes offset, monthly thin plate spline stratified by State, year*state interaction, random effect of IBGE6
# this function also checks fitting time
s = Sys.time()
m_base <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re"),
              family = "nb",
              method = "REML",
              data = ccc)
e = Sys.time()


# on my computer this takes 10 seconds - much better!
e-s 

################################################

# test lags for temperature
# adding 1 climate variable at a time to the baseline model

m_t0 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean, bs="tp"),
              family = "nb",
              method = "REML",
              data = ccc)

m_t1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_t2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_t3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_t4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_t5 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_t6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(tmean_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

# check AICs
aic_check = data.frame(
  lag = 0:6,
  AIC = unlist(lapply(list(m_t0, m_t1, m_t2, m_t3, m_t4, m_t5, m_t6), AIC))
)

#-------------------------------PRECIPITATION----------------------------------
m_p0 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p5 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_p6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + s(IBGE6, bs="re") + s(pmonth_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

#-------------------------------SPI-6----------------------------------
m_s0 <- gam(NumCases ~ -1 + offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp")+ s(spi6, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s1 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_1m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s2 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_2m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s3 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_3m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s4 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_4m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s5 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_5m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

m_s6 <- gam(NumCases ~ offset(offset) + s(Month, bs="tp", by=State) + Year*State + 
              s(IBGE6, bs="re") + s(pmonth, bs="tp")+ s(tmean_1m, bs="tp") + s(spi6_6m, bs="tp"),
            family = "nb",
            method = "REML",
            data = ccc)

AIC(m_base, m_t0, m_t1, m_t2, m_t3, m_t4, m_t5, m_t6, 
    m_p0, m_p1, m_p2, m_p3, m_p4, m_p5, m_p6, m_s0, m_s1, m_s2, m_s3, m_s4, m_s5, m_s6)

# extract and plot spline from best fitting model
# n.b. apply exponential transform to gain incidence rate ratio
draw(m_s0, fun=exp, select="s(tmean_1m)", caption=FALSE)+ 
  labs(title="(A) Association between lagged Temperature and Disease Incidence Rate", 
       x="Monthly Mean Temperature with 1-month lag (°C)", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="chikungunya-T1.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s0, fun=exp, select="s(pmonth)", caption=FALSE)+ 
  labs(title="(B) Association between Precipitation and Disease Incidence Rate", 
       x="Monthly Precipitation (mm)", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  scale_x_continuous(limits=c(0,400))+
  scale_y_continuous(limits=c(0,6))+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="chikungunya-P_short.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s0, fun=exp, select="s(spi6)", caption=FALSE)+ 
  labs(title="(C) Association between SPI-6 and Disease Incidence Rate", 
       x="SPI-6", y="Incidence Rate Ratio") +
  geom_hline(yintercept=1, lty=2, color="grey40")+
  scale_x_continuous(limits=c(-6,2))+
  scale_y_continuous(limits=c(0,3))+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename="chikungunya-S_short.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s0, fun=exp, select=c("s(Month):StateMato Grosso") , caption=FALSE) + 
  labs(title="Seasonal Variation of Chikungunya Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="MG.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s0, fun=exp, select=c("s(Month):StateGoiás") , caption=FALSE) + 
  labs(title="Seasonal Variation of Chikungunya Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="GO.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

draw(m_s0, fun=exp, select=c("s(Month):StateMato Grosso do Sul") , caption=FALSE) + 
  labs(title="Seasonal Variation of Chikungunya Incidence Rate",
       y="Incidence Rate Ratio") +
  scale_x_continuous(breaks=seq(1,12,1))+
  geom_hline(yintercept=1, lty=2, color="grey40")+
  theme_minimal()
ggsave(filename="MS.jpg", device="jpg", units="in", width=10, height=4, dpi=600)

#========================================================================================
# Example: assign rough population sizes (or use actual means from your data)
nn<-cc |> dplyr::select(IBGE6, State, Region, Municipality, Population, Year, Month)

# create state-year population table
state_pop <- nn %>%
  filter(Month==1) %>%
  group_by(State, Year) %>%
  summarise(State_Pop = sum(Population)) |> 
  ungroup()

print(state_pop)

new_data <- new_data |> filter(State%in%c("Mato Grosso do Sul", "Mato Grosso", "Goiás")) 

# Extract model predictions for year-month combinations
new_data <- expand.grid(
  Month = 1:12,
  Year = 2017:2024,
  State = c("Goiás", "Mato Grosso", "Mato Grosso do Sul"),
  IBGE6 = "dummy", # placeholder
  tmean_1m = mean(ccc$tmean_1m),
  pmonth = mean(ccc$pmonth),
  spi6 = mean(ccc$spi6)
)

# Merge with new_data by State and Year
new_data <- merge(new_data, state_pop, by = c("State", "Year"))

# Create the Offset
new_data$offset <- log(new_data$State_Pop / 100000)

# Predict and plot
predictions <- predict(m_s0, newdata=new_data, type="response")
new_data$predicted <- predictions

ggplot(new_data, aes(x=Month, y=Year, fill=predicted)) +
  geom_tile() +
  facet_wrap(~State) +
  scale_fill_viridis_c(name="Fitted\nIncidence rate") +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=seq(2017,2024,1)) +
  theme_minimal() +
  labs(title="Fitted chikungunya incidence by month and year across states")

# ------------------------------Diagnosis-------------------------------------
# check the output
summary(m_s0) 
 

# Extract year effects with CIs
year_coefs <- coef(m_s0)[grep("Year", names(coef(m_s0)))]
year_se <- summary(m_s0)$se[grep("Year", names(coef(m_s0)))]

# Calculate CIs
year_ci_lower <- year_coefs - 1.96 * year_se
year_ci_upper <- year_coefs + 1.96 * year_se

# Get first 8 coefficients
year_coefs_main <- year_coefs[1:8]
year_ci_lower_main <- year_ci_lower[1:8]
year_ci_upper_main <- year_ci_upper[1:8]

# Create results table
year_results <- data.frame(
  Year = 2017:2024,
  Log_IRR = year_coefs_main,
  CI_Lower = year_ci_lower_main,
  CI_Upper = year_ci_upper_main,
  IRR = exp(year_coefs_main),
  IRR_CI_Lower = exp(year_ci_lower_main),
  IRR_CI_Upper = exp(year_ci_upper_main)
)
print(year_results)

# check residuals
par(mfrow = c(2, 2))  # Set up 2 rows and 2 columns
gam.check(m_s4)
ggsave(filename="gam_check.jpg", device="jpg", units="in", width=5, height=4, dpi=600)


save(list = ls(all.names = TRUE), file = "chikungunya.RData") 
