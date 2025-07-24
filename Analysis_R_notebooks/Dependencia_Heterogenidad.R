# Proyecto Análisis Geoespacial
# Dependencia espacial (SAR y CAR)
# Base de datos de conteo de Balso (GBIF data) a nivel de Municipio

# Librerías SAR
library(sf)         # Para trabajar con datos espaciales (simple features)
library(spdep)      # Para matrices de vecindad y tests de dependencia espacial
library(spatialreg) # Para los modelos de regresión espacial
library(dplyr)      # Para manipulación de datos (similar a pandas)
# Manejo de base de datos
require(magrittr)
require(janitor)
require(tidyverse)
require(kbl)
require(kableExtra)

setwd("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data")
gdf = st_read("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data/Balso_Municipios.gpkg")

# Limpiar los nombres
gdf %<>% clean_names()
gdf %>% names

# Aplicar transformació log + 1
gdf <- gdf %>%
  mutate(y = log(numpoints + 1))

# Definir las variables dependiente e independientes
dependent_var <- "y"
independent_vars <- c("elev_mean",
                      "temperatura_media_anual_mean",
                      "precipitacion_anual_mean",
                      "rango_medio_diurno_mean",
                      "precipitacion_mes_mas_lluvioso_mean",
                      "precipitacion_mes_mas_seco_mean",
                      "isotermalidad_mean",
                      "estacionalidad_de_la_temperatura_mean",
                      "rango_anual_de_temperatura_mean",
                      "estacionalidad_de_la_precipitacion_mean")

# Crear la fórmula para el modelo
formula <- as.formula(paste(dependent_var, "~", paste(independent_vars, collapse = " + ")))

# Spatial matrix
W.nb <- poly2nb(gdf)
# Eliminar áreas sin vecinos
has_neighbors <- card(W.nb) > 0
gdf <- gdf[has_neighbors, ]
W.nb_clean <- subset(W.nb, subset = has_neighbors)
# Convertir a matriz y lista de pesos espaciales
W.mat <- nb2mat(W.nb_clean, style="B", zero.policy=TRUE)
W.list <- nb2listw(W.nb_clean, style="B", zero.policy=TRUE)
# Convertir para INLA (falta el archivo)
nb2INLA(file = "gdf.graph", nb = W.nb_clean)

# Cálculo de la matriz de vecindad (Spatial Weights Matrix, con k-vecinos más cercanos)
# Se obtienen los centroides de cada polígono/punto para calcular las distancias
# db_coords <- st_coordinates(st_centroid(gdf))
# Matriz de vecindad de k-vecinos más cercanos
# db_nb <- knn2nb(knearneigh(db_coords, k = 112)) #k=número de vecinos más cercanos.
# Convertir la matriz de vecindad a un objeto de tipo ‘listw’ (list of weights)
# db_listw <- nb2listw(db_nb, style = "W")

#-------------------------------------------------------------------------------
# MODELOS SAR
#-------------------------------------------------------------------------------
# Esta metodología de Dependencia espacial se ajustan a un regresión lm

"Notas: Las funciones lagsarlm(), lmSLX() y errorsarlm() asumen errores normales y
no está diseñado para familias GLM (como binomial)"

#-------------------------------------------------------------------------------
# 1. Spatially Lagged X (SLX)
#-------------------------------------------------------------------------------
slx_model <- lmSLX(formula, data = gdf, listw = W.list)
summary(slx_model)

# Moran's I for residuals of SAR model
slx_residuals <- residuals(slx_model)
moran_slx <- moran.mc(slx_residuals, listw = W.list, nsim = 999, zero.policy = TRUE)
print(moran_slx)

#-------------------------------------------------------------------------------
# 2. Spatial Autoregressive Model (SAR)
#-------------------------------------------------------------------------------
sar_model <- lagsarlm(formula, data = gdf, listw = W.list, type = "lag", zero.policy = TRUE)
summary(sar_model)

# Effects for SAR model
sar_effects <- impacts(sar_model, listw = W.list, R = 200, zero.policy = TRUE, trace = FALSE)
summary(sar_effects, zstats = TRUE, short = TRUE)

# Moran's I for residuals of SAR model
sar_residuals <- residuals(sar_model)
moran_sar <- moran.mc(sar_residuals, listw = W.list, nsim = 999, zero.policy = TRUE)
print(moran_sar)

#-------------------------------------------------------------------------------
# 3. Spatial Error Model (SEM)
#-------------------------------------------------------------------------------
sem_model <- errorsarlm(formula, data = gdf, listw = W.list, zero.policy = TRUE)
summary(sem_model)

# Moran's I for residuals of SEM model
sem_residuals <- residuals(sem_model)
moran_sem <- moran.mc(sem_residuals, listw = W.list, nsim = 999, zero.policy = TRUE)
print(moran_sem)

#-------------------------------------------------------------------------------
# MODELOS SAR + HETEROGENIDAD
#-------------------------------------------------------------------------------
# apply the log + 1 transformation
gdf <- gdf %>%
  mutate(y_transf = log(conteo + 1))

# Regímenes - Modelo de Retardo Espacial (SAR) - cuenca
model_sar_heter <- lagsarlm(y_transf ~ 0 + (elev_mean + 
                                       Temperatura_media_anual_mean + 
                                       Precipitacion_anual_mean + 
                                       Rango_medio_diurno_mean +
                                       Precipitacion_mes_mas_lluvioso_mean + 
                                       Precipitacion_mes_mas_seco_mean +
                                       Isotermalidad_mean + 
                                       Estacionalidad_de_la_temperatura_mean +
                                       Rango_anual_de_temperatura_mean +
                                       Estacionalidad_de_la_precipitacion_mean):(departamento), 
                     data = gdf, listw = W.list)

summary(model_sar_heter, Nagelkerke=T)

#-------------------------------------------------------------------------------
# MODELOS CAR
#-------------------------------------------------------------------------------
# Librerías
library(spBayes)
library(maps)
library(RANN)
library(gjam)
library(CARBayes)
library(CARBayesdata)
library(mgcv)
library(spdep)

# Heterogenidad y dependencia
library(sf) # Simple features for spatial data
library(spdep) # Spatial dependence tools
library(spatialreg) # Spatial regression models
library(sf)
library(dplyr)
library(sp)
library(spdep)
library(INLA) 

# Analisis CAR con datos a nivel de Municipios
# gdf = st_read("G:/My Drive/Investigacion2025/Posgrado_Statistics/GeoAnalysis/data/Balso_Municipios.gpkg")
# gdf %<>% clean_names()

# Variables dependientes e independientes
# Crear la variable binaria 'evento'
gdf <- gdf %>% mutate(evento = if_else(numpoints >= 1, 1, 0))
dependent_var <- "evento"

independent_vars <- c("elev_mean",
                      "temperatura_media_anual_mean",
                      "precipitacion_anual_mean",
                      "rango_medio_diurno_mean",
                      "precipitacion_mes_mas_lluvioso_mean",
                      "precipitacion_mes_mas_seco_mean",
                      "isotermalidad_mean",
                      "estacionalidad_de_la_temperatura_mean",
                      "rango_anual_de_temperatura_mean",
                      "estacionalidad_de_la_precipitacion_mean")

gdf <- gdf %>% mutate_at(independent_vars, ~(scale(.) %>% as.vector))

# Crear la fórmula para el modelo
formula <- as.formula(paste(dependent_var, "~", paste(independent_vars, collapse = " + ")))

#-------------------------------------------------------------------------------
# Standard Binomial model in GLM,  SIN DEPENDENCIA ESPACIAL 
#-------------------------------------------------------------------------------
# Ajustar modelo logístico clásico
model_glm <- glm(formula, data = gdf, family = binomial(link = "logit"))

# Resumen del modelo
summary(model_glm)

#-------------------------------------------------------------------------------
# Standard Binomial model in CARBayes, SIN DEPENDENCIA ESPACIAL 
#-------------------------------------------------------------------------------
model_carbayes <- S.glm(formula, data=gdf, family="binomial", trials = rep(1, nrow(gdf)),  
                     burnin=10000, n.sample=30000, thin=10,
                     n.chains=2, n.cores=1)

saveRDS(model_carbayes, file = "model_carbayes.rds") #Guardar el caché del modelo

model_carbayes

model_carbayes$modelfit %>%
  kbl(digits = 3, caption = "Resumen del modelo binomial S.glm") %>%
  kable_styling(full_width = FALSE)

model_carbayes$summary.results %>%
  kbl(digits = 3, caption = "Resumen del modelo binomial S.glm") %>%
  kable_styling(full_width = FALSE)

model_carbayes <- readRDS("model_carbayes.rds")

#-------------------------------------------------------------------------------
# Standard Binomial model in ICAR, CON DEPENDENCIA ESPACIAL (rho=1)
#-------------------------------------------------------------------------------
model_car_icar <- S.CARleroux(formula, data=gdf, family="binomial", trials=rep(1, nrow(gdf)),
                        W=W.mat, rho=1, 
                        burnin=10000, n.sample=30000, thin=10,
                        n.chains=2, n.cores=1)

saveRDS(model_car_icar, file = "model_car_icar.rds") #Guardar el caché del modelo

model_car_icar
model_car_icar$modelfit

#-------------------------------------------------------------------------------
# Standard Binomial model in ICAR, CON DEPENDENCIA ESPACIAL (BYM)
#-------------------------------------------------------------------------------
model_car_bym <- S.CARbym(formula, data=gdf, family="binomial", trials=rep(1, nrow(gdf)),
                   W=W.mat,
                   burnin=10000, n.sample=30000, thin=10,
                   n.chains=2, n.cores=1)
saveRDS(model_car_bym, file = "model_car_bym.rds") #Guardar el caché del modelo

model_car_bym
model_car_bym$modelfit

#-------------------------------------------------------------------------------
# ¿Cómo comparar S.glm(), S.CARleroux() y S.CARbym() ? Usa los criterios bayesianos de ajuste:
#-------------------------------------------------------------------------------
# 1. DIC (Deviance Information Criterion)
# DIC más bajo = mejor ajuste, penaliza la complejidad del modelo
# 2. p.d (número de parámetros efectivos)
# Da una idea de la complejidad del modelo
# 3. LMPL (Log Marginal Pseudo-Likelihood)
# Más alto (menos negativo) = mejor ajuste
model_carbayes$modelfit
model_car_icar$modelfit
model_car_bym$modelfit

#-------------------------------------------------------------------------------
# INLA. 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Modelo INLA Regresión logística SIN DEPENDENCIA ESPACIAL
#-------------------------------------------------------------------------------
modelo_inla <- inla(formula, data = as.data.frame(gdf),family = "binomial",
                control.predictor = list(compute = TRUE),
                control.compute = list(dic = TRUE, waic = TRUE))
summary(modelo_inla)
saveRDS(modelo_inla, file = "modelo_inla.rds") #Guardar el caché del modelo

#-------------------------------------------------------------------------------
# Modelo INLA Regresión logística CON HETEROGENEIDAD ESPACIAL
#-------------------------------------------------------------------------------
modelo_inla_heter <- inla(
  formula = evento ~ elev_mean + temperatura_media_anual_mean + precipitacion_anual_mean + 
    rango_medio_diurno_mean + precipitacion_mes_mas_lluvioso_mean + 
    precipitacion_mes_mas_seco_mean + isotermalidad_mean + estacionalidad_de_la_temperatura_mean + 
    rango_anual_de_temperatura_mean + estacionalidad_de_la_precipitacion_mean + 
    f(dpto_cnmbr, model="iid"),
  data = as.data.frame(gdf),
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(modelo_inla_heter)
saveRDS(modelo_inla_heter, file = "modelo_inla_heter.rds") #Guardar el caché del modelo

#-------------------------------------------------------------------------------
# Modelo INLA Regresión logística basado en modelo (besag)
#-------------------------------------------------------------------------------
nb2INLA(file = "gdf.graph", nb = W.nb_clean)

gdf$id <- 1:nrow(gdf)  # Asegúrate de que esté bien indexado

## CON DEPENDENCIA ESPACIAL---------------------------------------------------
modelo_inla_besag <- inla(
  formula = evento ~ elev_mean + temperatura_media_anual_mean + precipitacion_anual_mean + 
    rango_medio_diurno_mean + precipitacion_mes_mas_lluvioso_mean + 
    precipitacion_mes_mas_seco_mean + isotermalidad_mean + estacionalidad_de_la_temperatura_mean + 
    rango_anual_de_temperatura_mean + estacionalidad_de_la_precipitacion_mean + 
    f(id, model = "besag", graph = "gdf.graph"),
  data = as.data.frame(gdf),
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(modelo_inla_besag)
saveRDS(modelo_inla_besag, file = "modelo_inla_besag.rds") #Guardar el caché del modelo

## CON DEPENDENCIA Y HETEROGENEIDAD ESPACIAL ---------------------------------
modelo_inla_besag_heter_depend <- inla(
  formula = evento ~ elev_mean + temperatura_media_anual_mean + precipitacion_anual_mean + 
    rango_medio_diurno_mean + precipitacion_mes_mas_lluvioso_mean + 
    precipitacion_mes_mas_seco_mean + isotermalidad_mean + estacionalidad_de_la_temperatura_mean + 
    rango_anual_de_temperatura_mean + estacionalidad_de_la_precipitacion_mean + 
    f(dpto_cnmbr, model="iid")+
    f(id, model = "besag", graph = "gdf.graph"),
  data = as.data.frame(gdf),
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(modelo_inla_besag_heter_depend)
saveRDS(modelo_inla_besag_heter_depend, file = "modelo_inla_besag_heter_depend.rds") #Guardar el caché del modelo

## Comparar entre modelos-----------------------------------------------------
summary(modelo_inla)
summary(modelo_inla_heter)
summary(modelo_inla_besag)
summary(modelo_inla_besag_heter_depend)

## Mapa de probabilidades de ocurrencia del evento 
library(sf)
library(ggplot2)
library(viridis)  # Para una paleta de colores amigable

# Añadir las probabilidades de predicción al GeoDataFrame
gdf$probability_evento <- modelo_inla_besag_heter_depend$summary.fitted.values$mean

# Mapa de predicción del mejor modelo
inla_besag <- ggplot(gdf) +
  geom_sf(aes(fill = probability_evento), color = "grey80", size = 0.2) +
  scale_fill_viridis(name = "Probabilidad de Evento", option = "C", limits = c(0, 1)) +
  labs(title = "Mapa de Probabilidad de Ocurrencia del Evento") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )


#-------------------------------------------------------------------------------
# Modelo INLA Regresión logística CON DEPENDENCIA ESPACIAL (BYM)
#-------------------------------------------------------------------------------
modelo_inla_bym <- inla(
  formula = evento ~ elev_mean + temperatura_media_anual_mean + precipitacion_anual_mean + 
    rango_medio_diurno_mean + precipitacion_mes_mas_lluvioso_mean + 
    precipitacion_mes_mas_seco_mean + isotermalidad_mean + estacionalidad_de_la_temperatura_mean + 
    rango_anual_de_temperatura_mean + estacionalidad_de_la_precipitacion_mean + 
    f(id, model = "bym", graph = "gdf.graph"),
  data = as.data.frame(gdf),
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(modelo_inla_bym)
saveRDS(modelo_inla_bym, file = "modelo_inla_bymg.rds") #Guardar el caché del modelo

## CON DEPENDENCIA Y HETEROGENEIDAD ESPACIAL ---------------------------------
modelo_inla_bym_heter_depend <- inla(
  formula = evento ~ elev_mean + temperatura_media_anual_mean + precipitacion_anual_mean + 
    rango_medio_diurno_mean + precipitacion_mes_mas_lluvioso_mean + 
    precipitacion_mes_mas_seco_mean + isotermalidad_mean + estacionalidad_de_la_temperatura_mean + 
    rango_anual_de_temperatura_mean + estacionalidad_de_la_precipitacion_mean + 
    f(dpto_cnmbr, model="iid")+
    f(id, model = "bym", graph = "gdf.graph"),
  data = as.data.frame(gdf),
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(modelo_inla_bym_heter_depend)
saveRDS(modelo_inla_bym_heter_depend, file = "modelo_inla_bym_heter_depend.rds") #Guardar el caché del modelo

## Comparar entre modelos-----------------------------------------------------
model_carbayes <- readRDS("model_carbayes.rds")
model_car_icar <- readRDS("model_car_icar.rds")
model_car_bym <- readRDS("model_car_bym.rds")

modelo_inla <- readRDS("modelo_inla.rds")
modelo_inla_heter <- readRDS("modelo_inla_heter.rds")
modelo_inla_besag <- readRDS("modelo_inla_besag.rds")
modelo_inla_besag_heter_depend <- readRDS("modelo_inla_besag_heter_depend.rds")
modelo_inla_bym <- readRDS("modelo_inla_bym.rds")

summary(model_glm) #Modelo logístico clásico. Sin dependencia ni heterogenidad espacial
model_carbayes
model_car_icar
model_car_bym
summary(modelo_inla) #Modelo logístico usando INLA
summary(modelo_inla_heter)
summary(modelo_inla_besag)
summary(modelo_inla_besag_heter_depend)
summary(modelo_inla_bym)
summary(modelo_inla_bym_heter_depend)

## Mapa de probabilidades de ocurrencia del evento 
library(sf)
library(ggplot2)
library(viridis)  # Para una paleta de colores amigable

# Añadir las probabilidades de predicción al GeoDataFrame
gdf$probability_evento <- modelo_inla_besag_heter_depend$summary.fitted.values$mean

# Mapa de predicción del mejor modelo
inla_besag <- ggplot(gdf) +
  geom_sf(aes(fill = probability_evento), color = "grey80", size = 0.2) +
  scale_fill_viridis(name = "Probabilidad de Evento", option = "C", limits = c(0, 1)) +
  labs(title = "Mapa de Probabilidad de Ocurrencia del Evento") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )


# Añadir las probabilidades de predicción al GeoDataFrame
gdf$probability_evento_bym <- modelo_inla_bym_heter_depend$summary.fitted.values$mean

inla_bym <- ggplot(gdf) +
  geom_sf(aes(fill = probability_evento_bym), color = "grey80", size = 0.2) +
  scale_fill_viridis(name = "Probabilidad de Evento", option = "C", limits = c(0, 1)) +
  labs(
    title = "Mapa de Probabilidad de Ocurrencia del Evento",
    subtitle = "Predicciones generadas a partir del modelo INLA BYM + Heterogeneidad"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )



