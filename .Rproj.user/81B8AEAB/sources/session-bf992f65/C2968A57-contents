#install.packages("pheatmap")

#dependencias
library(readxl)
library(dplyr)
library(naniar)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(FSA)
library(dunn.test)

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
mrna

colnames(datos_geneticos) %>% 
  tail(n = 175) #la base de datos dice que hay 175 pero parece haber 173

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
colnames(datos_mrna) # selecciono los mrna (z scores)

#genes asociados a cáncer de mama: BRCA1, BRCA2, PALB2, CHEK2, CDH1
# Asume que tu data frame de Z-Scores se llama 'datos_expresion'
genes_buscados <- c("brca1", "brca2", "palb2", "chek2", "cdh1")

# Usar grep para buscar esos nombres en las columnas
genes_encontrados_expresion <- colnames(datos_mrna) %>%
  grep(paste(genes_buscados, collapse = "|"), ., value = TRUE, ignore.case = TRUE)

print("Genes encontrados en la tabla de Expresión (Z-Scores):")
print(genes_encontrados_expresion) #(5/5)
#esos genes están en mi tabla, vamos a ver si están relacionados con cáncer


# y a ver si están en mutaciones
genes_encontrados_expresion_mut <- colnames(datos_mutaciones) %>%
  grep(paste(genes_buscados, collapse = "|"), ., value = TRUE, ignore.case = TRUE)

print("Genes encontrados en la tabla de Expresión (Z-Scores):")
print(genes_encontrados_expresion_mut)
#hay mutaciones de ch1, brca1, brca2 y chek2

##############################3
# GRÁFICOS #################
#heatmap
pheatmap(datos_mrna[, genes_encontrados_expresion])

#relación entre z score brca2 y mutation_count
grafico_mutationz <- data %>% 
  select(brca2, mutation_count) %>% 
  filter(!is.na(brca2) & !is.na(mutation_count))

# creamos el Gráfico de Dispersión
dispersion_brca2 <- grafico_mutationz %>%
  ggplot(aes(x = brca2, y = mutation_count)) +
  # Puntos de dispersión (pacientes)
  geom_point(alpha = 0.6, color = "darkblue") +
  # Línea de regresión lineal (método lm) con su intervalo de confianza
  # Esto visualiza la correlación
  geom_smooth(method = "lm", color = "red", se = TRUE) + 
  labs(
    title = "Correlación entre Expresión de BRCA2 y Carga Mutacional",
    x = "Expresión de BRCA2 (Z-Score)",
    y = "Conteo de Mutaciones (Inestabilidad Genómica)"
  ) +
  theme_bw()

print(dispersion_brca2) # no parece haber relación

cor.test(grafico_mutationz$brca2, grafico_mutationz$mutation_count, method= "pearson")
# 0.03 corr

# grafico conteo de mutaciones frente subtipo
# Asumiendo que 'datos_combinados' está unido y limpio de NAs en las columnas
unique(data$pam50_._claudin.low_subtype)


datos_combinados <- data %>% 
  select(`pam50_._claudin.low_subtype`, mutation_count) %>% 
  filter(!is.na(`pam50_._claudin.low_subtype`) & !is.na(mutation_count)) 
  
  datos_combinados %>% ggplot(aes(x = pam50_._claudin.low_subtype, y = mutation_count, fill = pam50_._claudin.low_subtype)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) + # Violin plot muestra la densidad de los datos
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, show.legend = FALSE) + # Añade la mediana/quartiles
  labs(
    title = "Carga Mutacional por Subtipo Molecular PAM50",
    x = "Subtipo Molecular",
    y = "Conteo de Mutaciones",
    fill = "Subtipo" 
  ) +
  theme_bw() +
  scale_fill_viridis_d()
  
  
# vamos a hacerlo otra vez omitiendo NC
  # Limpieza y preparación de datos (Ajustado)
datos_combinados <- data %>% 
select(`pam50_._claudin.low_subtype`, mutation_count) %>% 
    
    # AÑADIR FILTRO PARA ELIMINAR "NC" (Not Classified)
filter(!is.na(`pam50_._claudin.low_subtype`) & 
             !is.na(mutation_count) & 
             `pam50_._claudin.low_subtype` != "NC") 
  
  # Nota: Después de filtrar, R puede conservar los niveles 'NC' en la metadata.
  # Para asegurarnos de que desaparezca del gráfico, convertimos la columna a factor, 
  # forzando a R a recalcular los niveles (se hace automáticamente si la columna ya era factor).
  datos_combinados <- datos_combinados %>%
    mutate(`pam50_._claudin.low_subtype` = factor(`pam50_._claudin.low_subtype`))
  
  
  # Gráfico (Box-Violin Plot)
  datos_combinados %>% 
    ggplot(aes(x = `pam50_._claudin.low_subtype`, y = mutation_count, fill = `pam50_._claudin.low_subtype`)) +
    geom_violin(alpha = 0.6, show.legend = FALSE) + 
    geom_boxplot(width = 0.1, color = "black", alpha = 0.7, show.legend = FALSE) + 
    labs(
      title = "Carga Mutacional por Subtipo Molecular PAM50",
      x = "Subtipo Molecular",
      y = "Conteo de Mutaciones",
      fill = "Subtipo" 
    ) +
    theme_bw() +
    scale_fill_viridis_d()

#existen diferencias significativas?
#voy a usar kruskal wallis
kruskal.test(mutation_count ~ `pam50_._claudin.low_subtype`, data= datos_combinados) 
#existen diferencias significativas

#para ver entre qué grupos existe diferencia: post hoc, Dunns
dunn_subtipo <- dunn.test(
  x= datos_combinados$mutation_count,
  g= datos_combinados$`pam50_._claudin.low_subtype`,
  method= "bonferroni",
  altp= TRUE, # muestra p valor ajustado
  table = TRUE,# muestra resultados en formato tabla
)


#filtramos solo datos significativos
