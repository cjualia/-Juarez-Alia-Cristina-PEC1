# Cargo SummarizedExperiment

load("SummarizedExperiment_GCdata.Rda")

dim(se)

# INSTALACIÓN DE PAQUETES NECESARIOS
# Biobase (SummarizedExperiment)
# ggplot2
# plotly
# scatterplot3d
# factoextra

# -----x-------x-------x-----PREPROCESAMIENTO------x---------x------------x-----

#-------------------------- LIMPIEZA DE DATOS---------------------------------

### Eliminamos del estudio aquellos metabolitos con QC-RSD <25%
### rowData() y colData() son independientes del assay. 
### Tenemos que modificar el assay directamente. 

se <-se[(rowData(se)$QC_RSD < 25) & (rowData(se)$Perc_missing < 20), ] 

dim(se) # Nos quedamos con 72 metabolitos

#------------------ TRATAMIENTO DE NA: IMPUTACIÓN POR VALOR MEDIO--------------

## Creamos una función que sustituya los NA por la media correspondiente a ese 
# metabolito, según el grupo

imputar_medias <- function(se){
  
  BN <- rownames(colData(se)[colData(se)$clase == "BN", ])
  HE <- rownames(colData(se)[colData(se)$clase == "HE", ])
  GC <- rownames(colData(se)[colData(se)$clase == "GC", ])
  QC <- rownames(colData(se)[colData(se)$clase == "QC", ])
  
  indices <- which(is.na(assay(se)), arr.ind = TRUE)
  
  metabolitos <- rownames(indices)
  
  ind_col <- indices[,"col"]
  muestras <- colnames(assay(se))[ind_col]
  
  for (i in seq_along(metabolitos)) {
    muestra <- muestras[i]
    metabolito <- metabolitos[i]
    
    clase_muestra <- colData(se)[muestra,]$clase
    
    if (clase_muestra == "BN"){
      media = mean(assay(se)[metabolito,BN], na.rm = TRUE)
    } else if (clase_muestra == "HE"){
      media = mean(assay(se)[metabolito,HE], na.rm = TRUE)
    } else if (clase_muestra == "GC"){
      media = mean(assay(se)[metabolito,GC], na.rm = TRUE)
    } else if (clase_muestra == "QC"){
      media = mean(assay(se)[metabolito,QC], na.rm = TRUE)
    }
    assay(se)[metabolito, muestra] = media
  }
  return(se)
}

## Aplicamos la función al SummarizedExperiment
se <- imputar_medias(se)

## Comprobamos que ya no existen NA:
anyNA(assay(se))

#---------------ANÁLISIS DE LOS QC: EVALUACIÓN DE LA REPRODUCIBILIDAD-----------

## Con los datos SIN ESCALAR (para no distorsionar el significado de los QCs)

pca1 <- prcomp(t(assay(se)), center = TRUE, scale. = TRUE)

library(dplyr)

datos_pca1 <- as.data.frame(pca1$x[, 1:3]) %>% 
  cbind(clase = colData(se)$clase, batch = colData(se)$batch)

pca_2d <- ggplot(datos_pca1) +
  aes(PC1, PC2, color = clase) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c('black', 'black', 'black', "red")) +  
  ggtitle("Análisis de QC") +
  theme(plot.title = element_text(hjust = 0.5))  

pca_2d

## Una vez hemos utilizado los QC para evaluar la calidad del experimento en
## términos de reproducibilidad; sustituimos el SummarizedEsperiment actual
## por otro sin datos de QC

QC <- rownames(colData(se)[colData(se)$clase == "QC", ])
nueva_matrix <- assay(se)[, !(colnames(assay(se)) %in% QC), drop = FALSE]
nuevo_colData <- as.data.frame(colData(se)[!(colnames(assay(se)) %in% QC), , drop = FALSE])
nuevo_colData$clase <- droplevels(nuevo_colData$clase)

se <- SummarizedExperiment(
  assays = SimpleList(counts = nueva_matrix), 
  colData = nuevo_colData,                      
  rowData = rowData(se)                      
)

dim(se) # Ahora tiene 123 muestras

rm(QC, nueva_matrix, nuevo_colData)

####### A LO MEJOR NO ES NECESARIO #####

#---------------------------ESCALAR DATOS---------------------------------------

## Escalamos los datos por el método del escalado Z
## Utilizaremos la media y desviación estándar correspondiente a cada grupo
## No escalamos los datos de los QC
## Diseñamos una función que devuelve una matriz con los datos escalados de 
## todos los metabolitos del grupo (clase) especificado

scale_by_group <- function(se, clase) {
  
  
  BN <- rownames(colData(se)[colData(se)$clase == "BN", ])
  HE <- rownames(colData(se)[colData(se)$clase == "HE", ])
  GC <- rownames(colData(se)[colData(se)$clase == "GC", ])
  
  
  if (clase == "BN") {
    muestras = BN
  } else if (clase == "HE") {
    muestras = HE
  } else if (clase == "GC") {
    muestras = GC
  } 
  
  # Calcular medias y desviaciones estándar
  medias = apply(assay(se)[, muestras], 1, mean, na.rm = TRUE)
  sds = apply(assay(se)[, muestras], 1, sd, na.rm = TRUE)
  
  # Matriz con datos escalados 
  m <- (assay(se)[,muestras] - medias)
  
  return(m) 
}

# Utilizamos lapply para utilizar como argumento 'clase' cada una de las 
# denominaciones de cada grupo (BN, HE, GC, QC)
# El resultado es una lista de las matrices resultantes de la función
results <- lapply(c("BN","HE","GC"), function(clase) {
  scale_by_group(se, clase)
})


# Cada uno de los elementos de 'results' se pasa como argumento a cbind().
# Así conseguimos unir todas las matrices que contienen los datos normalizados
assay_scaled <- do.call(cbind, results)

# Sustituimos el assay(se) original por 'assay_scale'
# Primero reorganizamos el assay escalado
assay_scaled <- assay_scaled[, match(colnames(assay(se)), colnames(assay_scaled)), drop = FALSE]
assay(se) = assay_scaled

rm(assay_scaled, results)

# ---------------------------ANÁLISIS MULTIVARIANTE ----------------------------

## ---------------------------------PCA----------------------------------------

### Evaluamos varias cosas:
### Análisis de calidad:
### > Existe efecto batch? 
### > Reproducibilidad? Muestras QC se agrupan entre sí? SÍ!
### Análisis de agrupamiento
### > Las muestras se agrupan según clase? Es percibible mayor variación en muestras GC y BN, menos en HE
### > Metabolitos que más influyen en la separación de los grupos? 

# Realizar PCA sobre la matriz transpuesta de datos de assay(se)
# Indicamos scale = FALSE porque ya hemos escalado anteriormente
pca2 <- prcomp(t(assay(se)), center = TRUE, scale. = FALSE)

library(dplyr)

datos_pca2 <- as.data.frame(pca2$x[, 1:3]) %>% 
  cbind(clase = colData(se)$clase, batch = colData(se)$batch)

# Calculamos las varianzas explicadas por cada componente principal

var_expl1 <- round(((pca2$sdev[1])^2/sum((pca2$sdev)^2)*100),2)
var_expl2 <- round(((pca2$sdev[2])^2/sum((pca2$sdev)^2)*100),2)
var_expl3 <- round(((pca2$sdev[3])^2/sum((pca2$sdev)^2)*100),2)
var_expl_total <- sum(var_expl1,var_expl2,var_expl3)

# Crear el gráfico 3D con scatterplot3d

### Genero una función que me permita realizar dos gráficos PCA:
### > Para evaluar el efecto batch
### > Para evaluar la variabilidad entre grupos

library(scatterplot3d)

plot_pca <- function(df, colores, factor, titulo) {
  # Crear el gráfico 3D
  scatterplot3d(df$PC1, df$PC2, df$PC3, 
                color = colores[as.numeric(factor)],  
                pch = 20,  
                xlab = paste("PC1", var_expl1, "%"),  
                ylab = paste("PC2", var_expl2, "%"),  
                zlab = paste("PC3", var_expl3, "%"))
  
  # Añadir la leyenda
  legend("topright",                             
         legend = levels(factor),   
         col = colores[1:length(levels(factor))],  
         pch = 19,
         inset = c(0.001, 0.001)) 
  
  # Añadir el título
  title(main = titulo)
}

### Aplico la función:
par(mfrow = c(1,2))

#### PCA según batch

colores2 <- c("#0D3B66", "#1D5E8A", "#377BA7", "#A3D5E0")

plot_pca(df=datos_pca2, colores=colores2, factor=datos_pca2$batch,
         titulo=paste("PCA según batch"))
### PCA según clase

colores1 <- c('deeppink', 'darkviolet', 'green3')

plot_pca(df=datos_pca2, colores=colores1, factor=datos_pca2$clase,
         titulo=paste("PCA según grupo de muestra"))

mtext("Varianza explicada total:", var_expl_total, "%")

par(mfrow = c(1,1))

### Evaluamos los 10 metabolitos más influyentes en la distribución de variabilidad entre grupos:

metabo_PC1 <- sort(abs(pca2$rotation[,1]), decreasing = TRUE)
metabo_PC1_r <- metabo_PC1[1:15]
nombres <- rowData(se)[names(metabo_PC1_r), ]$Label


## BOXPLOT MULTIVARIANTE

# grouped boxplot

## Diseño una función que me devuelva los valores de una serie de metabolitos de una clase

valores_met <- function(se,metabolitos,clase){
  
  BN <- rownames(colData(se)[colData(se)$clase == "BN", ])
  HE <- rownames(colData(se)[colData(se)$clase == "HE", ])
  GC <- rownames(colData(se)[colData(se)$clase == "GC", ])
  
  if (clase == "BN"){
    muestras = BN
  } else if (clase == "HE"){
    muestras = HE
  } else if (clase == "GC"){
    muestras = GC
  }
  
  df = data.frame()
  for (metabolito in metabolitos){
   datos = data.frame(assay(se)[metabolito, muestras])
   datos = cbind(datos,rep(clase,length(muestras)), rep(metabolito,length(muestras)))
   df = rbind(df, datos)
  }
  colnames(df) = c("valor", "clase", "metabolito")
  return(df)
}


# Aplica `valores_met()` para cada clase 
seleccion <- lapply(c("BN", "HE", "GC"), function(clase)
  {valores_met(se = se, metabolitos = names(metabo_PC1_r), clase = clase)})

# Combina los resultados en un solo dataframe,lo utilizaremos para el boxplot
df_boxplot <- do.call(rbind, seleccion)

# grouped boxplot
ggplot(df_boxplot, aes(x=metabolito, y=valor, fill=clase)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("BN" = "deeppink", "GC" = "darkviolet", "HE" = "green3"))+
  ggtitle("Metabolitos más influyentes en PC1")


##-----------------------------CLUSTERING--------------------------------------

# Seleccionamos solo dos grupos (por simplificar el análisis)

GC <- rownames(colData(se)[colData(se)$clase == "GC", ])
HE <- rownames(colData(se)[colData(se)$clase == "HE", ])

# Creamos la nueva matriz con los datos seleccionados

matrix = assay(se)[names(metabo_PC1),c(HE,GC)]

# Vamos a añadir al nombre de cada muestra la terminación BN, GC o HE
# en función del grupo al que pertenezca.
# Esto facilita la interpretación de los próximos gráficos.

nombres_columnas <- colnames(matrix)
clases_muestras <- colData(se)[nombres_columnas, ]$clase
nuevos_nombres <- paste0(nombres_columnas, "_", clases_muestras)
colnames(matrix) = substr(nuevos_nombres, 8, nchar(nuevos_nombres))

### DENDOGRAMA: AGRUPAMIENTO DE MUESTRAS

dist_matrix <- dist(t(matrix), method = "euclidean")

hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc, main = "Dendrograma de Clustering Jerárquico", xlab = "Muestras", cex = 0.7)

### HEATMAP DE MUESTRAS (otra representación de lo anterior)

heatmap(
  matrix,
  Rowv = TRUE,
  Colv = TRUE,
  distfun = function(x) dist(x, method = "euclidean"),
  hclustfun = function(x) hclust(x, method = "ward.D2"),
  scale = "row",
  col = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap con Clustering de Muestras",
  cexRow = 0.5,  # Ajusta el tamaño de la fuente para las filas (muestras)
  cexCol = 0.5   # Ajusta el tamaño de la fuente para las columnas (variables)
)

### MÉTODO DE K-MEANS

library(factoextra)

fviz_nbclust(t(matrix), kmeans, method = "wss") + 
  labs(title = "")

km <- kmeans(t(matrix), centers = 2, iter.max = 1000)


fviz_cluster(km, data = t(matrix), 
             geom = "point", 
             stand = FALSE, 
             ellipse.type = "convex", 
             main = "K-means Clustering") +
  geom_text(aes(label = rownames(t(matrix))), 
            size = 3, 
            vjust = 1, 
            hjust = 1) + 
  theme_minimal()





