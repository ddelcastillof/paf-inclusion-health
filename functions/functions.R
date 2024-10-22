#########################
###### FUNCTIONS ########
#########################

# To calculate attributable deaths and population attributable fractions

attr_deaths <- function(prop_inc_health, smr, deaths) 
{
  (prop_inc_health * (smr - 1)) / (prop_inc_health * (smr - 1) + 1) * deaths
}


paf <- function(prop_inc_health, smr) 
{
  (prop_inc_health * (smr - 1)) / (prop_inc_health * (smr - 1) + 1)
}