#instalar paquetes necesarios
#install.packages("dplyr")
#install.packages("naniar")
#install.packages("FSA")
#install.packages("survival")
#install.packages("survminer")
#install.packages("splines")

#dependencias
library(readxl)
library(dplyr)
library(naniar)
library(ggplot2)
library(tidyr)
library(FSA)
library(survival)
library(survminer)
library(splines)

#cargamos base de datos
data <- read.csv("C:/Users/usuario/Downloads/archive/METABRIC_RNA_Mutation.csv")
View(data)

###### 1. Exploracion inicial de datos #########################
dim(data) #1904 filas 693 columnas

colnames(data) %>% 
  head(n = 40) #vemos las primeras 40 columnas

colnames(data) %>% 
  tail(n = 40) #las 40 últimas

summary(data)

glimpse(data)

anyDuplicated(colnames(data)) #para ver si hay columnas replicadas


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


#vamos a ver los datos perdidos
# Detectar valores vacíos
sapply(datos_clinicos, function(x) sum(x == "", na.rm= TRUE))

# converir valores vacíos a NA
datos_clinicos[datos_clinicos == ""] <- NA 

#NA por columna
na_clinicos <- colSums(is.na(datos_clinicos)) #te devuelve la tabla entera
na_clinicos[na_clinicos>1] #para ver qué variables tienen más NA

# Revisar porcentaje de NA por columna
colMeans(is.na(datos_clinicos)) * 100
which(colMeans(is.na(datos_clinicos))
       > 0.02) #tumor_stage es la columna con mayor número de datos perdidos, con el 30% de datos perdidos, en total 501.


# Porcentaje de NA por fila ??
porc_na_fila <- rowMeans(is.na(datos_clinicos)) * 100
porc_na_fila #asi no

#resumen estadístico rápido
summary(datos_clinicos)

#para manejar los datos faltantes??
#eliminamos filas con NA, o usamos media para imputar o no los tenemos en cuenta a la hora de hacer test/modelos

# Quiero hacer gráficos de "tumor_stage"
#creo un df para hacer el gráfico
str(datos_clinicos$tumor_stage) #es numérico
datos_clinicos$tumor_stage <- factor(datos_clinicos$tumor_stage) #hacemos que sea de tipo categórico
boxplot(datos_clinicos$tumor_stage)

boxplot_supervivencia <- datos_clinicos %>%
  # X es la etapa, Y es la supervivencia en meses
  ggplot(aes(x = tumor_stage, y = overall_survival_months, fill = tumor_stage)) +
  geom_boxplot() +
  # Añado puntos para ver la distribución individual
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  labs(
    title = "Supervivencia Global por Etapa del Tumor (METABRIC)",
    x = "Etapa del Tumor",
    y = "Supervivencia Global (Meses)",
    fill = "Etapa" # Leyenda del relleno
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis_d(option= "D") # Usamos 'fill' para el relleno de las cajas, colourblind friendly

print(boxplot_supervivencia)

#Ahora quitamos los NA
datos_grafico <- datos_clinicos %>%
  filter(!is.na(tumor_stage)) %>% #filtramos NA 
  mutate(tumor_stage = factor(tumor_stage)) #hacemos que sea tipo factor

# 3. Verificar si hay NAs en la columna antes de graficar
sum(is.na(datos_grafico$tumor_stage)) #da 0

boxplot_supervivenciafinal <- datos_grafico %>%
  # X es la etapa, Y es la supervivencia en meses
  ggplot(aes(x = tumor_stage, y = overall_survival_months, fill = tumor_stage)) +
  geom_boxplot(show.legend= FALSE) +
  # Añado puntos para ver la distribución individual
  geom_jitter(color="black", size=0.4, alpha=0.5, show.legend = FALSE) +
  labs(
    title = "Supervivencia Global por Etapa del Tumor (METABRIC)",
    x = "Etapa del Tumor",
    y = "Supervivencia Global (Meses)",
    fill = "Etapa" # Leyenda del relleno
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis_d(option= "D")

print(boxplot_supervivenciafinal)

# Existen diferencias significativas?

#Comprobamos normalidad
shapiro.test(datos_grafico$overall_survival_months) # distribución no normal
kruskal.test(overall_survival_months ~ tumor_stage, data= datos_grafico) #existen diferencias

#hacemos post-hoc para ver entre qué etapas hay diferencia
resultado_dunn <- dunnTest(overall_survival_months ~ tumor_stage,
         data= datos_grafico, method= "bonferroni")

tabla_dunn <- resultado_dunn$res #extraemos tabla de resultados
significativos <- tabla_dunn %>% 
  filter(P.adj < 0.05) #seleccionamos solo los significativamente diferentes
significativos %>% 
  select(Comparison, P.adj, Z)

#la tasa de supervivencia es mayor en la etapa 1, es significativamente diferente de la 2, 3 y 4.
#diferencia de etapa 2 a 3, mejor tasa de supervivencia en la 2.
#A mayor etapa de tumor, peor supervivencia.

#quiero graficar la edad de diagnóstico y la supervivencia

#quiero graficar el tamaño del tumor y la supervivencia
coxph(Surv(overall_survival_months, overall_survival) ~ tumor_size, data = datos_clinicos) #no hay relación entre tamaño del tumor y supervivencia


ggsurvplot(survfit(Surv(overall_survival_months, overall_survival) ~ tumor_stage, 
                   data = datos_clinicos)) #parece que hay diferencias

ggplot(datos_clinicos, aes(x = tumor_size, y = overall_survival_months)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Relación entre tamaño del tumor y supervivencia",
       x = "Tamaño del tumor (mm)",
       y = "Supervivencia global (meses)") +
  theme_minimal()

#otra prueba


