#ANALISIS DISCRIMINANTE

library(MVN)
library(biotools)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(MASS)
library(pROC)
library(caTools)

datos <- read.csv("datos_limpios.csv")

#---------------NORMALIDAD------------------------

vars <- c("slc_pfp_track_score",
          "slc_pfp_shower_energy",
          "slc_pfp_track_chi2_muon",
          "slc_pfp_shower_length",
          "slc_nu_score",
          "slc_pfp_shower_opening_angle")

clases <- unique(datos$target)

for(clase in clases){
  
  cat("\n=====================\n")
  cat("TARGET =", clase, "\n")
  cat("=====================\n")
  
  datos_clase <- datos[datos$target == clase, vars]
  
  n <- min(nrow(datos_clase), 1999)
  muestra <- datos_clase[sample(nrow(datos_clase), n), ]
  
  # Shapiro-Wilk por variable
  cat("\n--- Shapiro-Wilk ---\n")
  for(v in vars){
    cat("\nVariable:", v, "\n")
    print(shapiro.test(muestra[[v]]))
  }
  
  # Henze-Zirkler
  cat("\n--- Henze-Zirkler ---\n")
  hz <- mvn(data = muestra, mvnTest = "hz")
  print(hz$multivariateNormality)
  
  # Royston
  cat("\n--- Royston ---\n")
  roy <- mvn(data = muestra, mvnTest = "royston")
  print(roy$multivariateNormality)

}  



#---------Q-Q plots por variable y clase ---------------------
for(clase in clases){
  datos_clase <- datos[datos$target == clase, vars]
  muestra <- datos_clase[sample(nrow(datos_clase), 500), ]
  
  par(mfrow = c(2, 3))
  for(v in vars){
    qqnorm(muestra[[v]], main = paste(v, "-", clase))
    qqline(muestra[[v]], col = "red")
  }
}


#-------------------HOMOCEDASTICIDAD------------------------

# Test M de Box
box_test <- boxM(datos[, vars], datos$target)
print(box_test)

for(clase in clases){
  cat("\nMatriz de covarianzas - TARGET =", clase, "\n")
  muestra <- datos[datos$target == clase, vars]
  print(round(cov(muestra), 3))
}


#-----------------------APLICAMOS EL MODELO---------------------------

#Dividimos los datos en 70% entrenamiento y 30% test
split <- sample.split(datos$target, SplitRatio = 0.7)

train_data <- subset(datos, split == TRUE)
test_data  <- subset(datos, split == FALSE)

#LDA 
modelo_lda <- lda(target ~ slc_pfp_track_score + slc_pfp_shower_energy +
                    slc_pfp_track_chi2_muon + slc_nu_score + slc_pfp_shower_length +
                    slc_pfp_shower_opening_angle,
                  data = train_data) #conjunto de entrenamiento

# Predicciones
probs <- predict(modelo_lda, test_data)$posterior[, 2]
pred_clase <- ifelse(probs > 0.5, 1, 0)

# Matriz de confusion
matriz_conf <- table(Real = test_data$target, Predicho = pred_clase)
print(matriz_conf)

VN <- matriz_conf[1, 1]
FP <- matriz_conf[1, 2]
FN <- matriz_conf[2, 1]
VP <- matriz_conf[2, 2]

# Metricas
eficiencia  <- VP / (VP + FN)
especificidad <- VN / (VN + FP)
pureza        <- VP / (VP + FP)
producto      <- eficiencia * pureza
accuracy      <- (VP + VN) / (VP + VN + FP + FN)



# Curva ROC
roc_obj <- roc(test_data$target, probs)

roc_data <- data.frame(
  TPR = roc_obj$sensitivities,
  FPR = 1 - roc_obj$specificities
)

grafico <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#008B8B", linewidth = 1) +
  geom_area(fill = "#AFEEEE", alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.65, y = 0.35,
           label = paste("AUC == ", round(auc(roc_obj), 4)),
           parse = TRUE, size = 6, color = "#008B8B", fontface = "bold") +
  labs(title = "ROC Curve - LDA",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray92")
  )

ggsave("roc_lda.png", plot = grafico, width = 8, height = 6, dpi = 300)

cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n")
cat("AUC =", auc(roc_obj), "\n")





#------------------FRONTERA LINEAL---------------------

# Coeficientes de la funcion discriminante
print(modelo_lda$scaling)

# Punto de corte
intercepto <- -0.5 * t(modelo_lda$scaling) %*% 
  (modelo_lda$means[1,] + modelo_lda$means[2,])
cat("Punto de corte:", intercepto, "\n")






#----------------------QDA------------------------

# Ajuste del modelo
modelo_qda <- qda(target ~ ., data = train_data[, c(vars, "target")])

# Predicciones
probs_qda <- predict(modelo_qda, test_data)$posterior[, 2]
pred_clase_qda <- ifelse(probs_qda > 0.5, 1, 0)

# Matriz de confusion
matriz_conf_qda <- table(Real = test_data$target, Predicho = pred_clase_qda)
print(matriz_conf_qda)

VN <- matriz_conf_qda[1, 1]
FP <- matriz_conf_qda[1, 2]
FN <- matriz_conf_qda[2, 1]
VP <- matriz_conf_qda[2, 2]

# Metricas
eficiencia  <- VP / (VP + FN)
especificidad <- VN / (VN + FP)
pureza        <- VP / (VP + FP)
producto      <- eficiencia * pureza
accuracy      <- (VP + VN) / (VP + VN + FP + FN)



# Curva ROC
roc_obj_qda <- roc(test_data$target, probs_qda)

roc_data_qda <- data.frame(
  TPR = roc_obj_qda$sensitivities,
  FPR = 1 - roc_obj_qda$specificities
)

grafico_qda <- ggplot(roc_data_qda, aes(x = FPR, y = TPR)) +
  geom_line(color = "#008B8B", linewidth = 1) +
  geom_area(fill = "#AFEEEE", alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.65, y = 0.35,
           label = paste("AUC == ", round(auc(roc_obj_qda), 4)),
           parse = TRUE, size = 6, color = "#008B8B", fontface = "bold") +
  labs(title = "ROC Curve - QDA",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray92")
  )

ggsave("roc_qda.png", plot = grafico_qda, width = 8, height = 6, dpi = 300)

cat("AUC =", auc(roc_obj_qda), "\n")
cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n")


# Medias por clase
modelo_qda$means

# Matrices de covarianzas por clase
modelo_qda$scaling
