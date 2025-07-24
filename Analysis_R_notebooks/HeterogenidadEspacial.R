# Modelo SAR y heterogenidad espacial

# Manejo de base de datos
require(magrittr)
require(janitor)
require(tidyverse)
require(kbl)
require(kableExtra)

# Librerías SAR
library(sf)         # Para trabajar con datos espaciales (simple features)
library(spdep)      # Para matrices de vecindad y tests de dependencia espacial
library(spatialreg) # Para los modelos de regresión espacial
library(dplyr)      # Para manipulación de datos (similar a pandas)

# Librerías Heterogenidad
library(sf) #para importar datos geoespaciales
library(lme4) # para el modelo
library(pscl) # para calcular los R2
library(MuMIn) #para calcular los R2 en el modelo multinivel 
library(ggplot2) # Plotting library
library(ggspatial) # For adding north arrow and scale to maps
library(sjPlot) #para graficar los effectos
library(dplyr) # for data manipulation
library(pROC) # for the ROC curve
library(broom)

#Cargar datos
setwd("C:/Users/Valentina Cardona/GitHub/GeoAnalysis_2025/data")
data= st_read("C:/Users/Valentina Cardona/GitHub/GeoAnalysis_2025/data/Balso_Municipios_Macrocuenca.gpkg")

# Renombramos algunas variables
data <- data %>% rename(
  departamento = dpto_cnmbr,
  municipio = mpio_cnmbr,
  area = mpio_narea,
  conteo = NUMPOINTS
)

# Create Y_bin as 1 if lands_rec is 1 or greater, otherwise 0
data$Y_bin <- ifelse(data$conteo >= 1, 1, 0)

# Convert 'departamento' to a factor (categorical variable)
data$departamento <- as.factor(data$departamento)
data$Macrocuenca <- as.factor(data$Macrocuenca)

# Ensure predictor variables are standardized
data$elev_mean_std <- scale(data$elev_mean)
data$Temp_media_std <- scale(data$Temperatura_media_anual_mean)
data$Prec_media_std <- scale(data$Precipitacion_anual_mean)
data$Rango_medio_std <- scale(data$Rango_medio_diurno_mean)
data$Precipitacion_mes_lluvioso_std <- scale(data$Precipitacion_mes_mas_lluvioso_mean)
data$Precipitacion_mes_seco_std <- scale(data$Precipitacion_mes_mas_seco_mean)
data$Isotermalidad_std <- scale(data$Isotermalidad_mean)
data$Estacionalidad_temp_std <- scale(data$Estacionalidad_de_la_temperatura_mean)
data$Rango_anual_temp_std <- scale(data$Rango_anual_de_temperatura_mean)
data$Estacionalidad_prec_std <- scale(data$Estacionalidad_de_la_precipitacion_mean)

# MODELO DE REGRESIÓN LOGÍSTICO SIN EFECTOS ALEATORIOS ------------------------
# Fit the logistic regression model without random effects

m1 <- glm(
  Y_bin ~ elev_mean_std + Temp_media_std + Prec_media_std + Rango_medio_std +
    Precipitacion_mes_seco_std + Isotermalidad_std + Estacionalidad_temp_std +
    Estacionalidad_prec_std,
  data = data, 
  family = binomial(link = "logit")
)

# Print the summary of the model
summary(m1)

# Calculate pseudo-R^2 values
pseudo_r2 <- pR2(m1)
pseudo_r2

# Fitted values for model m3
data$m1fitted <- fitted(m1) 

# Extract model coefficients as a tidy data frame
model_summary <- tidy(m1) %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  rename(
    Variables = term,
    Estimate = estimate,
    Std_Error = std.error,
    Z_value = statistic,
    P_value = p.value
  )

# MODELO DE REGRESIÓN LOGÍSTICA CON INTERCEPTO VARIABLE (EFECTO ALEATORIO)-----
# El modelo m2 es un modelo multinivel o modelo jerárquico que reconoce la estructura 
# anidada de los datos, donde las observaciones individuales (municipios) están 
# agrupadas dentro de unidades más grandes (macrocuenca). 

# DEPARTAMENTO--------------------------
m2 <- glmer(
  Y_bin ~ elev_mean_std + Temp_media_std + Prec_media_std + Rango_medio_std +
    Precipitacion_mes_seco_std + Isotermalidad_std + Estacionalidad_temp_std +
    Estacionalidad_prec_std +
    (1 | departamento),
  data = data,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# Print the summary of the model
summary(m2)

# Fitted values
data$m2 <- fitted(m2)

# Extract random intercepts
ranef_intercept <- ranef(m2)$departamento


# MACROCUENCA--------------------------
m3 <- glm(
  Y_bin ~ elev_mean_std + Temp_media_std + Prec_media_std + Rango_medio_std +
    Precipitacion_mes_seco_std + Isotermalidad_std + Estacionalidad_temp_std +
    Estacionalidad_prec_std + Macrocuenca,
  data = data,
  family = binomial(link = "logit"))

# Print the summary of the model
summary(m3)
