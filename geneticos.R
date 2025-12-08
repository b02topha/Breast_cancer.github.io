#dependencias
library(readxl)
library(dplyr)
library(naniar)
library(ggplot2)
library(tidyr)
library(FSA)

#cargamos base de datos
data <- read.csv("C:/Users/usuario/Downloads/archive/METABRIC_RNA_Mutation.csv")

#######################################
# 2.2 DATOS GENÉTICOS: ARNm y MUTACIONES
###################################
# Datos genéticos: ARNm + mutaciones
datos_geneticos <- data %>% 
  select(32:693) #sabemos que el resto son genéticos

colgeneticos <- colnames(datos_geneticos) %>% 
  head(n= 693)
colgeneticos

mrna <- colnames(datos_geneticos) %>% 
  head(n= 520)
zscore

mutaciones <- colnames(datos_geneticos) %>% 
  tail(n = 173) #columnas de mutaciones
mutaciones 

#convertimos espacios vacíos en NAs
datos_geneticos[datos_geneticos == ""] <- NA

#número de NA por columna
na_geneticos <- colSums(is.na(datos_geneticos))
na_geneticos[na_geneticos > 0] # no hay NAs

#quiero extraer las columnas de mutaciones
# Seleccionar columnas que terminan en "_mut"
mutaciones_cols <- grep("_mut$", colnames(datos_geneticos), value = TRUE) #son 173

# Crear un dataframe solo con mutaciones
datos_mutaciones <- datos_geneticos %>% select(all_of(mutaciones_cols))

datos_mrna <- data %>% 
  select(32:520)
colnames(mrna) # selecciono los mrna (z scores)

#genes asociados a cáncer de mama: BRCA1, BRCA2, PALB2, CHEK2, CDH1
# Asume que tu data frame de Z-Scores se llama 'datos_expresion'
genes_buscados <- c("brca1", "brca2", "palb2", "chek2", "cdh1")

# Usar grep para buscar esos nombres en las columnas
genes_encontrados_expresion <- colnames(datos_expresion) %>%
  grep(paste(genes_buscados, collapse = "|"), ., value = TRUE, ignore.case = TRUE)

print("Genes encontrados en la tabla de Expresión (Z-Scores):")
print(genes_encontrados_expresion)
