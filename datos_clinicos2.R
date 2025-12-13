library(survival)
library(readxl)
library(dplyr)
library(survminer)

# vamos a explorar los datos con uniq

#cargamos base de datos
data <- read.csv("C:/Users/usuario/Downloads/archive/METABRIC_RNA_Mutation.csv")
View(data)

colnames(data) %>% 
  head(n = 40) #las 40 primeras

############## 2. Separación de datos #####################
#separamos columnas por tipo de dato 

# 2.1 DATOS CLÍNICOS ################################
#################################################
col_datosclinicos <- colnames(data) %>% 
  head(n = 31)
View(col_datosclinicos)
#sabemos por la información proporcionada con el dataset que hay datos clínicos y datos genéticos. Seleccionamos los datos clínicos

datos_clinicos <- data %>% 
  select(col_datosclinicos) #selecciono datos clinicos de toda la base de datos (data)

str(datos_clinicos)
dim(datos_clinicos)

unique(datos_clinicos$tumor_stage) # así puedo ver todos los valores de la columna, ver si hay NA.. etc
unique(datos_clinicos$tumor_other_histologic_subtype) #aquí vemos que no hay espacios en blanco

# quiero relacionar survival months con her2_status: el cáncer her2+ es más agresivo frente al her2-, mi hipótesis es que va a tener menor supervivencia los her2+
unique(datos_clinicos$her2_status)
str(datos_clinicos$her2_status)

#importante que no tenga NA y que sea un factor
datos_her2 <- datos_clinicos %>%
  filter(!is.na(her2_status)) %>%
  # Convertimos a factor, aunque es probable que ya lo sea
  mutate(her2_status = factor(her2_status))
str(datos_her2) #ahora si es un factor


unique(datos_clinicos$overall_survival)
str(datos_clinicos$overall_survival_months)

# 2. Crear el Objeto de Supervivencia (Surv Object)
# 'overall_survival_months' = tiempo
# 'overall_survival' = evento (1=fallecido, 0=vivo/censurado)
surv_object_her2 <- Surv(time = datos_her2$overall_survival_months, 
                         event = datos_her2$overall_survival)


## 2. Crear el objeto de supervivencia
# Tiempo: overall_survival_months
# Evento: overall_survival (1=fallecido, 0=vivo/censurado)
surv_object_her2 <- Surv(time = datos_her2_supervivencia$overall_survival_months, 
                         event = datos_her2_supervivencia$overall_survival)

# 3. Ajustar el Modelo de Kaplan-Meier
# Ajustamos la supervivencia en función del estado HER2
fit_km_her2 <- survfit(surv_object_her2 ~ her2_status, data = datos_her2)

# 4. Generar el Gráfico y el Log-rank Test
km_her2_plot <- ggsurvplot(fit_km_her2, 
                           data = datos_her2,
                           pval = TRUE,             # Muestra el valor p de la prueba Log-rank
                           risk.table = TRUE,       # Muestra la tabla de riesgo
                           legend.title = "Estado HER2",
                           title = "Curvas de Supervivencia Global por Estado HER2",
                           palette = c("#E7B800", "#2E9FDF") # Colores para HER2- y HER2+
)

print(km_her2_plot)

#ya que existen diferencias significativas entre tener her2+ y her2- vamos a hacer el modelo de cox
# 1. Ajustar el Modelo de Cox (Univariante)
# Usamos el estado HER2 como único predictor de la supervivencia
cox_her2_model <- coxph(surv_object_her2 ~ her2_status, data = datos_her2)

# 2. Imprimir el resumen del modelo para obtener el Hazard Ratio (HR)
summary(cox_her2_model)

# interpretación: El valor más importante es el Hazard Ratio = 1.3391.Los pacientes con cáncer de mama HER2 Positivo tienen un 33.91% más de riesgo de experimentar el evento (muerte) en cualquier momento dado que los pacientes con cáncer HER2 Negativo.

#hacemos un modelo multivariante para ver si el riesgo asociado a her2 se ve afectado por otros factores como: etapa del tumor, edad, quimioterapia.
# 1. Seleccionar y limpiar solo las filas necesarias para el modelo multivariante
datos_multivariante <- datos_clinicos %>%
  select(overall_survival_months, overall_survival, # Variables de tiempo/evento
         her2_status, tumor_stage, age_at_diagnosis, chemotherapy) %>% # Covariables
  
  # Filtrar NAs en CUALQUIERA de las variables del modelo
  filter(!is.na(her2_status) & 
           !is.na(tumor_stage) & 
           !is.na(age_at_diagnosis) & 
           !is.na(chemotherapy)) %>%
  
  # Asegurar que las variables categóricas sean factores
  mutate(her2_status = factor(her2_status),
         tumor_stage = factor(tumor_stage),
         chemotherapy = factor(chemotherapy)) 

# 2. Crear el objeto de supervivencia (Surv Object)
surv_object_multi <- Surv(time = datos_multivariante$overall_survival_months, 
                          event = datos_multivariante$overall_survival)

# ejectumamos modelo cox multivariante
# 3. Ajustar el Modelo de Cox Multivariante
cox_multivariante_model <- coxph(
  surv_object_multi ~ her2_status + tumor_stage + age_at_diagnosis + chemotherapy,
  data = datos_multivariante
)

# 4. Imprimir el resumen del modelo
summary(cox_multivariante_model)
# el estado her2+, en análisis univariante es pronóstico, pero no es un predictor en la supervivencia si se hace un modelo con quimioterapia y estadio tumoral

#existe colinealidad?
#calculamos vif en el modelo de COX
library(car)

# 1. Crear un modelo lineal con los mismos predictores
# Solo se usa para calcular el VIF de las covariables.
lm_model <- lm(overall_survival_months ~ her2_status + tumor_stage + age_at_diagnosis + chemotherapy, 
               data = datos_multivariante)

# 2. Calcular el Factor de Inflación de la Varianza (VIF)
vif_results <- vif(lm_model)

print("Resultados del Factor de Inflación de la Varianza (VIF):")
print(vif_results)
#que no haya colinealidad indica que el modelo está bien

###############################################
# AHORA QUIERO RELACIONAR LATERALIDAD DEL TUMOR CON SUPERVIVENCIA #
# 1. Seleccionar y limpiar los datos para incluir Lateralidad
unique(datos_clinicos$primary_tumor_laterality) #hay espacios en blanco que afectaran, hay que considerarlos como NAs
# 1. Seleccionar y limpiar los datos para incluir Lateralidad

# Creamos un nuevo data frame de análisis, filtrando los valores vacíos ("")
datos_analisis_lateralidad <- datos_clinicos %>%
  select(overall_survival_months, overall_survival, 
         primary_tumor_laterality) %>% 
  
  # FILTRO CLAVE: Excluir las filas donde primary_tumor_laterality es una cadena vacía ""
  filter(primary_tumor_laterality != "") %>% 
  
  # Asegurar que las variables categóricas sean factores
  mutate(primary_tumor_laterality = factor(primary_tumor_laterality)) 


# 2. Crear el objeto de supervivencia (Surv Object)
survival_object <- Surv(datos_analisis_lateralidad$overall_survival_months, 
                        datos_analisis_lateralidad$overall_survival)

# 3. Ajustar el modelo Kaplan-Meier por lateralidad
km_laterality <- survfit(survival_object ~ primary_tumor_laterality, 
                         data = datos_analisis_lateralidad)

# 4. Realizar el Test de Log-Rank
log_rank_test <- survdiff(survival_object ~ primary_tumor_laterality, 
                          data = datos_analisis_lateralidad)

print("Resultados del Test de Log-Rank para Lateralidad:")
print(log_rank_test)

# 5. Generar el gráfico de Kaplan-Meier
kaplan_meier_plot <- ggsurvplot(
  km_laterality, 
  data = datos_analisis_lateralidad, 
  pval = TRUE,                       # Mostrar el p-valor del Log-Rank
  risk.table = TRUE,                 # Mostrar la tabla de riesgo debajo del gráfico
  legend.title = "Lateralidad del Tumor",
  title = "Curvas de Supervivencia Global por Lateralidad del Tumor",
  xlab = "Tiempo en meses",
  palette = c("purple", "orange") # Colores para Izquierdo/Derecho
)

# 6. Mostrar el gráfico
print(kaplan_meier_plot)
#existe diferencia significativa

# vamos a hacer un modelo de Cox
# 1. Ajustar el Modelo de Cox (Univariante)
# Usamos primary_tumor_laterality como único predictor
cox_laterality_model <- coxph(survival_object ~ primary_tumor_laterality, 
                              data = datos_analisis_lateralidad)

# 2. Imprimir el resumen del modelo para obtener el Hazard Ratio (HR)
summary(cox_laterality_model)
#interpretación:Los pacientes con cáncer de mama en el lado Derecho tienen un 17.42% más de riesgo de morir en cualquier momento dado en comparación con los pacientes con cáncer de mama en el lado Izquierdo.

##########################################################################
# AHORA COMPARAMOS TRIPLE NEGATIVO para ver si es más agresivo y tiene peor pronóstico: se define por la ausencia de tres receptores, ER (receptor de estrógeno), PR (receptor de progesterona) y receptor HER2-
unique(datos_clinicos$er_status)
unique(datos_clinicos$pr_status)
unique(datos_clinicos$her2_status)

# Asumimos que "Negative" significa ausencia de receptor.
# Usaremos el data frame original de datos clínicos.
datos_tnbc <- datos_clinicos %>%
  select(overall_survival_months, overall_survival, 
         er_status, pr_status, her2_status) %>%
  
  # 1. Crear la variable TNBC_status
  mutate(
    TNBC_status = case_when(
      er_status == "Negative" & pr_status == "Negative" & her2_status == "Negative" ~ "TNBC",
      TRUE ~ "Non_TNBC"
    )
  ) %>%
  
  # 2. Filtrar filas con NAs en los marcadores de supervivencia
  filter(!is.na(overall_survival_months) & !is.na(overall_survival)) %>%
  
  # 3. Convertir a factor para el análisis
  mutate(TNBC_status = factor(TNBC_status))

# 4. Crear el objeto de supervivencia
surv_object_tnbc <- Surv(datos_tnbc$overall_survival_months, 
                         datos_tnbc$overall_survival)
# A. Kaplan-Meier (Ver la diferencia visual y obtener el p-valor)
km_tnbc <- survfit(surv_object_tnbc ~ TNBC_status, data = datos_tnbc)

ggsurvplot(km_tnbc, 
           data = datos_tnbc,
           pval = TRUE,             # Mostrar p-valor
           risk.table = TRUE,       # Mostrar tabla de riesgo
           legend.title = "Estado Triple Negativo",
           title = "Supervivencia Global (TNBC vs. No-TNBC)")

# B. Modelo de Cox (Cuantificar el riesgo HR)
cox_tnbc_model <- coxph(surv_object_tnbc ~ TNBC_status, data = datos_tnbc)

print("Resumen del Modelo de Cox para TNBC:")
summary(cox_tnbc_model)

#interpretación: es significativo, "El cáncer de mama Triple Negativo (TNBC) se asocia con un riesgo de muerte 21.24% mayor en cualquier momento dado en comparación con otros subtipos de cáncer de mama (No-TNBC)

#####################################
# quiero ver cuantas pacientes son pre-post menopáusicas
# 1. Limpieza y preparación
# Crear un data frame limpio solo con la variable de interés, eliminando valores vacíos o NAs.
datos_grafico_menopausia <- datos_clinicos %>%
  filter(!is.na(inferred_menopausal_state) & inferred_menopausal_state != "") %>%
  mutate(inferred_menopausal_state = factor(inferred_menopausal_state))

# 2. Generar el gráfico de barras
grafico_barras_menopausia <- ggplot(datos_grafico_menopausia, 
                                    aes(x = inferred_menopausal_state, fill = inferred_menopausal_state)) +
  geom_bar() +
  
  # Añadir las etiquetas de conteo encima de cada barra
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 5) +
  
  labs(title = "Distribución de Pacientes por Estado Menopáusico",
       x = "Estado Menopáusico",
       y = "Número de Pacientes") +
  
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(grafico_barras_menopausia)
