library(sf)
library(magrittr)
library(corrplot)
library(tidyverse)

setwd("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data")
data <- st_read("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data/Balso_Municipios.gpkg")

data %>% as.data.frame() %>%
  select(4:14) %>%
  na.omit() %>% 
  cor() %>% 
  round(2) %>% 
  corrplot(method = "circle", 
           type = "lower", 
           addCoef.col = "black", 
           number.cex = 0.7,   # Tamaño de los coeficientes de correlación
           tl.cex = 0.5,       # Tamaño del texto de las etiquetas de las variables
           diag = FALSE)

