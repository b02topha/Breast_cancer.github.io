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
