## https://bit.ly/3kbMN4W
## Manipulating and analyzing data with dplyr
## Load the tidyverse packages
library(tidyverse)

## Import data
## Data on animal species diversity and weights found 
## within plots at our study site.
surveys <- read_csv("data_raw/portal_data_joined.csv")

## inspect the data
str(surveys)

## preview the data
View(surveys)

## Select variables
## Lets take a closer look at just species_id & weight 
select(surveys, plot_id, species_id, weight)

## Drop variables
select(surveys,-record_id,-species_id)

## Filter
## Lets look at data only from the year 1995
filter(surveys, year == 1995)

## Stringing data operations
## Retain species_id, sex & weight for 
## species with weight < 5
surveys2 <- filter(surveys, weight < 5)
surveys_sml <- select(surveys2, species_id, sex, weight)

surveys_sml <-
  select(
  filter(surveys, weight < 5), species_id, sex, weight
  )

## Pipes
surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)

surveys_sml <- surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)

surveys_sml

## Challenge
## Using pipes, subset the surveys data to include 
## animals collected before 1995 and retain only the 
## columns year, sex, and weight.

## Mutate
## To create a new column of weight in kg:
surveys %>%
  mutate(weight_kg = weight / 1000)

## Two variables
surveys %>%
  mutate(
    weight_kg = weight / 1000,
    weight_lb = weight_kg * 2.2
  )

surveys %>%
  mutate(weight_kg = weight / 1000) %>%
  head()

surveys %>%
  filter(!is.na(weight)) %>%
  mutate(weight_kg = weight / 1000) %>%
  head()

## Challenge
## Create a new data frame from the surveys data that 
## meets the following criteria: 

## Contains only the species_id column and a new column 
## called hindfoot_cm, containing the hindfoot_length 
## values converted to centimeters.

## In this hindfoot_cm column, there should be no NA's 
## and all values should be less than 3.

## Split-apply-combine data analysis
## The summarize function
surveys %>%
  group_by(sex) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE)) %>%
  tail()

surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight))

surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight)) %>%
  print(n = 15)

surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(
    mean_weight = mean(weight),
    min_weight = min(weight)
  )

surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(
    mean_weight = mean(weight),
    min_weight = min(weight)) %>%
  arrange(min_weight)

surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(
    mean_weight = mean(weight),
    min_weight = min(weight)) %>%
  arrange(desc(mean_weight))

## Counting
surveys %>%
  count(sex)

surveys %>%
  group_by(sex) %>%
  summarise(count = n())

surveys %>%
  count(sex, sort = TRUE)

surveys %>%
  count(sex, species)

surveys %>%
  count(sex, species) %>%
  arrange(species, desc(n))

## Challenge
## 1. How many animals were caught in each plot_type 
## surveyed?
## 2. Use group_by() and summarize() to find the mean, 
## min, and max hindfoot
## length for each species (using species_id). Also add 
## the number of observations (hint: see ?n).
## 3. What was the heaviest animal measured in each year? 
## Return the columns year, genus, species_id, and weight.


## Reshaping with pivot_longer & pivot_wider
## Pivot wider
surveys_gw <- surveys %>%
  filter(!is.na(weight)) %>%
  group_by(plot_id, genus) %>%
  summarize(mean_weight = mean(weight))

str(surveys_gw)

surveys_wide <- surveys_gw %>%
  pivot_wider(
    names_from = genus, 
    values_from = mean_weight
  )

str(surveys_spread)

## Replace missing values with 0
surveys_wide <- surveys_gw %>%
  pivot_wider(
    names_from = genus, 
    values_from = mean_weight,
    values_fill = 0
  )
  
## Pivot longer
surveys_long <- surveys_wide %>%
  pivot_longer(
    -plot_id,
    names_to = "genus",
    values_to = "mean_weight"
  )
  
str(surveys_long)

surveys_long <- surveys_wide %>%
  pivot_longer(
    Baiomys:Spermophilus,
    names_to = "genus",
    values_to = "mean_weight"
  )
  
## Exporting data
surveys_complete <- surveys %>%
  filter(
    !is.na(weight),           # remove missing weight
    !is.na(hindfoot_length),  # remove missing hindfoot_length
    !is.na(sex)               # remove missing sex
    )

## Extract the most common species_id
species_counts <- surveys_complete %>%
  count(species_id) %>% 
  filter(n >= 50)

## Only keep the most common species
surveys_complete <- surveys_complete %>%
  filter(species_id %in% species_counts$species_id)

## Export data
write_csv(surveys_complete, 
  path = "data/surveys_complete.csv")