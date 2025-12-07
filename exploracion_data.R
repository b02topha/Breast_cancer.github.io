#instalar paquetes necesarios
#install.packages("dplyr")

#dependencias
library(readxl)
library(dplyr)


#cargamos base de datos
data <- read.csv("C:/Users/usuario/Downloads/archive/METABRIC_RNA_Mutation.csv")
View(data)

#exploracion de datos
dim(data) #1904 filas 693 columnas
colnames(data) %>% 
  head(n = 40)

colnames(data) %>% 
  tail(n = 40)

colnames(data) %>% 
  tail(n = 175)

summary(data)

glimpse(data)

#separamos columnas por tipo de dato 
# datos clínicos
col_datosclinicos <- colnames(data) %>% 
  head(n = 31)
View(col_datosclinicos)
#sabemos por la información proporcionada con el dataset que hay datos clínicos y datos genéticos. Seleccionamos los datos clínicos

datos_clinicos <- data %>% 
  select(col_datosclinicos)
str(datos_clinicos)


datos_geneticos <- data %>% 
  select(32:693) #sabemos que el resto son genéticos
