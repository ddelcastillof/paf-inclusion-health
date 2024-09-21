# ----------
# PURPOSE:
#
# The aim of this file is to estimate global PAFs for inclusion health groups.
# This is work in progress with comments being improved and functions code being added continuously.
# 
# Data sources at present not fully referenced, but currently it reads in data from a google sheet on: 
# - inclusion health population sizes
# - SMRs for inclusion health groups based on a previous review https://doi.org/10.1016/S0140-6736(17)31869-X
# - other data required for the analyses e.g. GDB death estimates and population sizes
# 
# CONTRIBUTORS: Rob Aldridge
# 
# TO DOs: - fix missing homeless data for Argentina and Belgium.
#         - make the age groups for the SMRs more comparable
#         - make calendar years more comparable (separate years and decide what to do with the data, check all cause vs all cause mortality data)
#         - add better data on male and female numbers for each country as at the moment assumption is based on 
#           UK data and 80% M 20% F
#         - add in population numbers for each country by age and sex
#         - Model the inclusion health overlap using ? Poisson capture recapture methods and MCM simulation
#         - Model the inclusion health distribution by age and sex for each country
# ----------


# ----------
# set up libraries
# ----------
pacman::p_load(googlesheets4, metafor, tidyverse, skimr, mc2d, patchwork)


# ----------
# read in data
# ----------
googlesheets4::gs4_deauth()

input_smr_all <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=570599438", sheet = "input_mortality", col_types = "ccccccccdcddccdcdccccdddccccc??????")
input_country_region_codes <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_country_and_region_codes")
input_mortality_numbers <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_IHME-GBD_2019_DATA")
input_population_numbers <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_world_bank_population_totals")
input_country_name_clean <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_country_name_clean")
input_homeless_pop <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_homeless")
input_prison_pop <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_prison")
input_sud_pop <- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_sud")
input_hic_country_list<- read_sheet("https://docs.google.com/spreadsheets/d/1onGl6q4gkPI4a2dtUvEeqgQEee0NYceAHBBpyOOXs88/edit#gid=1722368264", sheet = "input_hic")

glimpse(input_country_region_codes)


#set the number of simulations
sim_run <- 1000

# disable scientific notation
options(scipen = 999)


# ----------
# clean the smr data
# ----------

#filter the data to only include non-excluded or duplicate SMR data
smr_filter_all_cause <- input_smr_all %>%
  filter(Outcome_measure=="SMR", ICD10_chapter_number==0, is.na(Duplicate), is.na(Exclude), !is.na(lci), !is.na(uci)) %>%
  select (lci, uci, Outcome, Population, Country, First_author, Year_of_publication, Female) %>%
  as.data.frame() %>%
  rename (sex = Female)

#create standard error estimates from confidence intervals from study
smr_filter_all_cause$sei <- (log(smr_filter_all_cause$uci/smr_filter_all_cause$lci))/3.92

# ----------
# clean the total mortality data
# ----------

glimpse(input_mortality_numbers)

high_income_countries_paf <- input_mortality_numbers %>%
  group_by(country_name, sex_id) %>%
  summarise(
    deaths = sum(val, na.rm = T),
    deaths_lci = sum(lower, na.rm = T),
    deaths_uci = sum(upper, na.rm = T),
    sim_run = sim_run
  )%>%
  filter (deaths >50000)%>%
  left_join(input_country_name_clean, by = "country_name")%>%
  filter (!is.na(alpha_3_country)) %>%
  ungroup() %>%
  mutate(
    sex = as.character (
      case_when ( sex_id == 1 ~ "Male",
                  sex_id == 2 ~ "Female",
                  sex_id == 3 ~ "Both"))
  ) %>%
  uncount(sim_run) %>%
  group_by(alpha_3_country, sex) %>%
  mutate(
    row_id = as.character(row_number()),
    id = str_c(alpha_3_country, row_id, sep = "-"),
    idsex = str_c(alpha_3_country, sex, row_id, sep = "-")
  ) %>%
  ungroup() %>%
  select (id, deaths, sex, idsex) 

# ----------
# clean the homeless population data
# ----------

glimpse(input_homeless_pop)

input_homeless_pop <- input_homeless_pop %>%
  select(alpha_3_country, homeless_pop, homeless_pop_lower, homeless_pop_upper) %>%
  mutate(sim_run = sim_run) %>%
  uncount(sim_run)  %>%
  rowwise() %>%
  mutate(
    homeless_mc_se = (homeless_pop_upper - homeless_pop_lower)/(2*1.96),
    homeless_lambda = rnorm(1, mean = homeless_pop, sd = homeless_mc_se),
    homeless_lambda = ifelse(homeless_lambda < 0, 0, homeless_lambda),
    homeless_mc_poisson = rpois(1, lambda = homeless_lambda)) %>%
  group_by(alpha_3_country) %>%
  mutate(
    row_id = as.character(row_number()),
    id = str_c(alpha_3_country, row_id, sep = "-")
  ) %>%
  ungroup() %>%
  select (id, homeless_mc_poisson) 


# ----------
# clean the prison population data
# ----------

glimpse(input_prison_pop)

input_prison_pop <- input_prison_pop  %>%
  left_join(input_country_name_clean, by = "country_name")  %>%
  right_join(input_hic_country_list, by = "alpha_3_country")   %>%
  filter (hic == 1) %>%
  select(alpha_3_country, prison_pop, prison_pop_lower, prison_pop_upper) %>% 
  mutate(sim_run = sim_run) %>%
  uncount(sim_run) %>%
  rowwise() %>%
  mutate(
    prison_mc_se = (prison_pop_upper - prison_pop_lower)/(2*1.96),
    prison_lambda = rnorm(1, mean = prison_pop, sd = prison_mc_se),
    prison_lambda = ifelse(prison_lambda < 0, 0, prison_lambda),
    prison_mc_poisson = rpois(1, lambda = prison_lambda)) %>%
  group_by(alpha_3_country) %>%
  mutate(
    row_id = as.character(row_number()),
    id = str_c(alpha_3_country, row_id, sep = "-")
  ) %>%
  ungroup() %>%
  select (id, prison_mc_poisson) 

# ----------
# clean the sud population data
# ----------

glimpse(input_sud_pop)

input_sud_pop <- input_sud_pop %>%
  group_by(country_name) %>%
  summarise(
    sud_pop = sum(val, na.rm = T),
    sud_pop_lower = sum(lower, na.rm = T),
    sud_pop_upper = sum(upper, na.rm = T),
  )%>%
  left_join(input_country_name_clean, by = "country_name")%>%
  filter (!is.na(alpha_3_country)) %>%
  left_join(input_hic_country_list, by = "alpha_3_country")   %>%
  filter (hic == 1) %>%
  select(alpha_3_country, sud_pop, sud_pop_lower, sud_pop_upper) %>% 
  mutate(sim_run = sim_run) %>%
  uncount(sim_run) %>%
  rowwise() %>%
  mutate(
    sud_mc_se = (sud_pop_upper - sud_pop_lower)/(2*1.96),
    sud_lambda = rnorm(1, mean = sud_pop, sd = sud_mc_se),
    sud_lambda = ifelse(sud_lambda < 0, 0, sud_lambda),
    sud_mc_poisson = rpois(1, lambda = sud_lambda)
  ) %>%
  group_by(alpha_3_country) %>%
  mutate(
    row_id = as.character(row_number()),
    id = str_c(alpha_3_country, row_id, sep = "-")
  ) %>%
  ungroup() %>%
  select (id, sud_mc_poisson) 


# ----------
# clean the total population data
# ----------

glimpse(input_population_numbers)

population_numbers <-  input_population_numbers %>% 
  pivot_longer(
    cols = `1960`:`2021`
  )%>% 
  rename(
    year = name,
    gen_pop_tot = value
  ) %>%
  filter(year == 2019)%>%
  mutate(sim_run = sim_run) %>%
  uncount(sim_run)%>%
  group_by(alpha_3_country) %>%
  mutate(
    row_id = as.character(row_number()),
    id = str_c(alpha_3_country, row_id, sep = "-")
  ) %>%
  ungroup() %>%
  select (id, gen_pop_tot) 

# ----------
# estimate total inclusion health for each country and calculate p
# ----------

inclusion_health_pop_est <- input_homeless_pop %>%
  left_join(input_prison_pop, by = "id") %>%
  left_join(input_sud_pop, by = "id") %>%
  left_join(population_numbers, by = "id") %>%
  rowwise() %>% 
  mutate(
    overlap_homeless = rpert(1, min=0.01,mode=.66,max=1),
    overlap_prison = rpert(1, min=0.01,mode=.66,max=1),
    inclusion_health_pop_tot = sum(homeless_pert*overlap_homeless, sud_pert, prison_pert*overlap_prison, na.rm = T),
    prop_inc_health = inclusion_health_pop_tot / gen_pop_tot
  )  

# ----------
# fit random-effects model to estimate the SMRs
# ----------

smr_filter_all_cause_both <- smr_filter_all_cause %>%
  filter(sex=="Both") 
  
res_smr_all_cause_both <- rma.uni(yi = log(smr_filter_all_cause_both$Outcome), sei=smr_filter_all_cause_both$sei, data = smr_filter_all_cause_both, method="REML")
log_smr.point_est_both <- as.numeric(exp(res_smr_all_cause_both$b))
log_smr.se_both <- as.numeric(res_smr_all_cause_both$se)


smr_filter_all_cause_male <- smr_filter_all_cause %>%
  filter(sex=="Male") 

res_smr_all_cause_male <- rma.uni(yi = log(smr_filter_all_cause_male$Outcome), sei=smr_filter_all_cause_male$sei, data = smr_filter_all_cause_male, method="REML")
log_smr.point_est_male <- as.numeric(exp(res_smr_all_cause_male$b))
log_smr.se_male <- as.numeric(res_smr_all_cause_male$se)


smr_filter_all_cause_female <- smr_filter_all_cause %>%
  filter(sex=="Female") 

res_smr_all_cause_female <- rma.uni(yi = log(smr_filter_all_cause_female$Outcome), sei=smr_filter_all_cause_female$sei, data = smr_filter_all_cause_female, method="REML")
log_smr.point_est_female <- as.numeric(exp(res_smr_all_cause_female$b))
log_smr.se_female <- as.numeric(res_smr_all_cause_female$se)


skim(high_income_countries_paf)

# ----------
# create the PAF dataset
# ----------

attr_deaths <- function(prop_inc_health, smr, deaths) 
  {
  (prop_inc_health * (smr - 1)) / (prop_inc_health * (smr - 1) + 1) * deaths
}

paf <- function(prop_inc_health, smr) 
{
  (prop_inc_health * (smr - 1)) / (prop_inc_health * (smr - 1) + 1)
}


high_income_countries_paf <- high_income_countries_paf %>%
  left_join(inclusion_health_pop_est, by = "id") %>%
  mutate(
    gen_pop_tot = ifelse(sex == "Male", gen_pop_tot /2, gen_pop_tot),
    gen_pop_tot = ifelse(sex == "Female", gen_pop_tot /2, gen_pop_tot),
    inclusion_health_pop_tot = ifelse(sex == "Male", inclusion_health_pop_tot * 0.8, inclusion_health_pop_tot),
    inclusion_health_pop_tot = ifelse(sex == "Female", inclusion_health_pop_tot * 0.2, inclusion_health_pop_tot),
    prop_inc_health = inclusion_health_pop_tot / gen_pop_tot
  ) %>%
    rowwise() %>% 
    mutate(
      smr = ifelse(sex == "Male", rnorm(1, mean=log_smr.point_est_male, sd = log_smr.se_male), 0),
      smr = ifelse(sex == "Female", rnorm(1, mean=log_smr.point_est_female, sd = log_smr.se_female), smr),
      smr = ifelse(sex == "Both", rnorm(1, mean=log_smr.point_est_both, sd = log_smr.se_both), smr),
      attr_deaths = attr_deaths(prop_inc_health, smr, deaths),
      paf = paf (prop_inc_health, smr),
      alpha_3_country = substr(id,1,3)
    ) %>% 
    ungroup() %>% 
    group_by(alpha_3_country, sex)%>% 
    summarise(
      paf_50 = quantile(paf, probs=(0.5)),
      paf_025 = quantile(paf, probs=(0.025)),
      paf_975 = quantile(paf, probs=(0.975)),
      prop_inc_health_50 = quantile(prop_inc_health, probs=(0.5)),
      prop_inc_health_025 = quantile(prop_inc_health, probs=(0.025)),
      prop_inc_health_975 = quantile(prop_inc_health, probs=(0.975)),
      smr_50 = quantile(smr, probs=(0.5)),
      smr_025 = quantile(smr, probs=(0.025)),
      smr_975 = quantile(smr, probs=(0.975)),
      attr_deaths_50 = quantile(attr_deaths, probs=(0.5)),
      attr_deaths_025 = quantile(attr_deaths, probs=(0.025)),
      attr_deaths_975 = quantile(attr_deaths, probs=(0.975)),
      deaths_50 = quantile(deaths, probs=(0.5))
    )%>% 
    mutate(
      not_attr_deaths_50 = deaths_50 - attr_deaths_50
    )


# PAF
  p1 <- high_income_countries_paf %>%
    filter(sex=="Male") %>%
    ggplot(aes(x=alpha_3_country, y=paf_50))+ 
    geom_bar( stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=alpha_3_country, ymin=paf_025, ymax=paf_975), width=0.4, colour="forestgreen", alpha=0.9, size=1.5) +
    ggtitle("Male")
  
  p2 <- high_income_countries_paf %>%
    filter(sex=="Female") %>%
    ggplot(aes(x=alpha_3_country, y=paf_50))+ 
    geom_bar( stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=alpha_3_country, ymin=paf_025, ymax=paf_975), width=0.4, colour="forestgreen", alpha=0.9, size=1.5) +
    ggtitle("Female")
  
  p3 <-  high_income_countries_paf %>%
    filter(sex=="Both") %>%
    ggplot(aes(x=alpha_3_country, y=paf_50))+ 
    geom_bar( stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=alpha_3_country, ymin=paf_025, ymax=paf_975), width=0.4, colour="forestgreen", alpha=0.9, size=1.5) +
    ggtitle("Both")

  
p1 + p2 + p3 + plot_layout(ncol = 3, nrow = 1)

  
# Proportion inclusion health
  high_income_countries_paf %>%
    filter(sex=="Both") %>%
    ggplot(aes(x=alpha_3_country, y=prop_inc_health_50))+ 
    geom_bar(  stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=alpha_3_country, ymin=prop_inc_health_025, ymax=prop_inc_health_975), width=0.4, colour="forestgreen", alpha=0.9, size=1.5) +
    ggtitle("Proportion inclusion health")
  
  
# SMR
  high_income_countries_paf %>%
    filter(sex=="Both") %>%
    ggplot( aes(x=alpha_3_country, y=smr_50))+ 
    geom_bar( stat="identity", fill="forestgreen", alpha=0.5) +
    geom_errorbar( aes(x=alpha_3_country, ymin=smr_025, ymax=smr_975), width=0.4, colour="forestgreen", alpha=0.9, size=1.5) +
    ggtitle("SMR")
  

#total and attributable deaths  
  p1 <- high_income_countries_paf %>%
    filter(sex=="Male") %>%
    select(alpha_3_country, not_attr_deaths_50, attr_deaths_50)%>% 
    pivot_longer(
      cols = not_attr_deaths_50:attr_deaths_50,
      names_to = "deaths",
      values_to = "counts"
    ) %>% 
    ggplot(aes(fill=deaths, y=alpha_3_country, x=counts)) + 
    geom_bar(position="stack", stat="identity")+
    ggtitle("Male")
  
  
  p2 <- high_income_countries_paf %>%
    filter(sex=="Female") %>%
    select(alpha_3_country, not_attr_deaths_50, attr_deaths_50)%>% 
    pivot_longer(
      cols = not_attr_deaths_50:attr_deaths_50,
      names_to = "deaths",
      values_to = "counts"
    ) %>% 
    ggplot(aes(fill=deaths, y=alpha_3_country, x=counts)) + 
    geom_bar(position="stack", stat="identity")+
    ggtitle("Female")
  
  
  p3 <- high_income_countries_paf %>%
    filter(sex=="Both") %>%
    select(alpha_3_country, not_attr_deaths_50, attr_deaths_50)%>% 
    pivot_longer(
      cols = not_attr_deaths_50:attr_deaths_50,
      names_to = "deaths",
      values_to = "counts"
    ) %>% 
    ggplot(aes(fill=deaths, y=alpha_3_country, x=counts)) + 
    geom_bar(position="stack", stat="identity")+
    ggtitle("Both")

p1 + p2 + p3 + plot_layout(ncol = 3)