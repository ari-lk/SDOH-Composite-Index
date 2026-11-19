library(readr)
library(dplyr)
library(stringr)

# SVI: keep only reliable & populated tracts
svi <- read_csv("California.csv", show_col_types = FALSE) %>%
  mutate(GEOID = as.character(FIPS)) %>%
  filter(E_TOTPOP > 0, F_TOTAL != -999) %>%   # only drop sentinel/unpopulated
  select(
    GEOID, F_TOTAL,
    E_TOTPOP, E_HU, E_HH,
    E_POV150, E_UNEMP, E_HBURD, E_NOHSDP, E_UNINSUR,
    E_AGE65, E_AGE17, E_DISABL, E_SNGPNT, E_LIMENG,
    E_MINRTY, E_MUNIT, E_MOBILE, E_CROWD, E_NOVEH, E_GROUPQ,
    RPL_THEMES, RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4
  )

cat("SVI tracts kept:", nrow(svi), "\n")

# PLACES: parse everything as character then numbers
places <- read_csv("PLACES Diabetes.csv", show_col_types = FALSE,
                   col_types = cols(.default = col_character())) %>%
  filter(StateAbbr == "CA", MeasureId == "DIABETES", DataValueTypeID == "CrdPrv") %>%
  transmute(
    GEOID        = LocationID,
    diabetes     = parse_number(Data_Value),
    diab_lo      = parse_number(Low_Confidence_Limit),
    diab_hi      = parse_number(High_Confidence_Limit),
    places_pop   = parse_number(TotalPopulation),
    places_pop18 = parse_number(TotalPop18plus),
    footnote     = Data_Value_Footnote_Symbol
  )

cat("PLACES rows:", nrow(places),
    "| unique tracts:", n_distinct(places$GEOID),
    "| counties:", n_distinct(substr(places$GEOID, 1, 5)), "\n")
cat("PLACES NA diabetes:", sum(is.na(places$diabetes)), "\n")

# merge
merged <- inner_join(svi, places, by = "GEOID")

cat("Matched:", nrow(merged), "\n")
cat("PLACES-only:", nrow(anti_join(places, svi, by = "GEOID")), "\n")
cat("F_TOTAL distribution in merged:\n"); print(table(merged$F_TOTAL))

# outcome audit
summary(merged$diabetes)
cat("NAs in merged diabetes:", sum(is.na(merged$diabetes)), "\n")
cat("Cor SVI vs diabetes:", round(cor(merged$RPL_THEMES, merged$diabetes, use = "complete.obs"), 3), "\n")

# Diagnostic: does dropping F_TOTAL>0 lose many tracts?
read_csv("California.csv", show_col_types = FALSE) %>%
  mutate(dropped = F_TOTAL > 0 | E_TOTPOP == 0) %>%
  count(dropped)
