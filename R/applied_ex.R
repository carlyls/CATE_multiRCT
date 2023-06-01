## Making plots for applied data example

library(tidyverse)

#age forest plot
ages <- data.frame(quart=c(1:4), estimate=c(1.536, 2.337, 2.388, 3.723),
                   SE=c(.993, .974, .991, .964))
ages %>%
  ggplot
