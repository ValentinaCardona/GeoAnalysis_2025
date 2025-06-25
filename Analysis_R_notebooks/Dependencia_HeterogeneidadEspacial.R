# Modelo de dependencia (SAR) y heterogenidad espacial

# Cargar librerías
library(spdep)
library(spatialreg)
library(sf) 
library(dplyr)
library(magrittr)

#Cargar datos
setwd("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data")
gdf= st_read("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data/Balso_Municipios.gpkg")

# Renombramos algunas variables
gdf <- gdf %>% rename(
  departamento = dpto_cnmbr,
  municipio = mpio_cnmbr,
  area = mpio_narea,
  conteo = NUMPOINTS
)

# apply the log + 1 transformation
gdf <- gdf %>%
  mutate(y_transf = log(conteo + 1))

# Define the dependent and independent variables
dependent_var <- "y_transf"
independent_vars <- c('elev_mean',
                      'Temperatura_media_anual_mean', 'Precipitacion_anual_mean',
                      'Rango_medio_diurno_mean', 'Precipitacion_mes_mas_lluvioso_mean',
                      'Precipitacion_mes_mas_seco_mean', 'Isotermalidad_mean',
                      'Estacionalidad_de_la_temperatura_mean',
                      'Rango_anual_de_temperatura_mean',
                      'Estacionalidad_de_la_precipitacion_mean')

# # Create the Queen contiguity weights matrix. Tenemos islas. Exploremos con vecinos más cercano
# queen_nb <- poly2nb(gdf) #Matriz queen. 
# queen_lw <- nb2listw(queen_nb, style = "W", zero.policy = TRUE)

# Create the (Spatial Weights Matrix) k-vecinos más cercanos
gdf_coords <- st_coordinates(st_centroid(gdf))
gdf_nb <- knn2nb(knearneigh(gdf_coords, k = 4))
gdf_listw <- nb2listw(gdf_nb, style = "W") #Convertir la matriz de vecindad a un objeto de tipo ‘listw’ (list of weights), que es el formato requerido


# Regímenes - Modelo de Retardo Espacial (SAR) - cuenca
col.fit9 <- lagsarlm(y_transf ~ 0 + (elev_mean + 
                                       Temperatura_media_anual_mean + 
                                       Precipitacion_anual_mean + 
                                       Rango_medio_diurno_mean +
                                       Precipitacion_mes_mas_lluvioso_mean + 
                                       Precipitacion_mes_mas_seco_mean +
                                       Isotermalidad_mean + 
                                       Estacionalidad_de_la_temperatura_mean +
                                       Rango_anual_de_temperatura_mean +
                                       Estacionalidad_de_la_precipitacion_mean):(departamento), 
                     data = gdf, listw = gdf_listw)

summary(col.fit9, Nagelkerke=T)

