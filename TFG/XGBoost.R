# XGBOOST


library(xgboost)
library(pROC)
library(ggplot2)
library(caTools)

datos<-read.csv("datos_limpios.csv")


set.seed(123) 
division <- sample.split(datos$target, SplitRatio = 0.70)

train_set <- subset(datos, division == TRUE)
test_set <- subset(datos, division == FALSE)



# XGBoost requiere matrices numéricas 
# Convertimos los sets a matrices
train_matrix <- as.matrix(train_set[, c("slc_pfp_track_score", "slc_pfp_shower_energy",
                                        "slc_pfp_track_chi2_muon", "slc_nu_score",
                                        "slc_pfp_shower_length", "slc_pfp_shower_opening_angle")])
test_matrix  <- as.matrix(test_set[, c("slc_pfp_track_score", "slc_pfp_shower_energy",
                                       "slc_pfp_track_chi2_muon", "slc_nu_score",
                                       "slc_pfp_shower_length", "slc_pfp_shower_opening_angle")])

# Entrenar el modelo 
modelo_xgb <- xgboost(data = train_matrix, 
                      label = as.numeric(as.character(train_set$target)), 
                      nrounds = 100, 
                      objective = "binary:logistic",
                      verbose = 0)

# Obtener probabilidades continuas
predicciones_prob <- predict(modelo_xgb, test_matrix)

# 1. Definir el umbral y convertir probabilidades a clases (0 o 1)
umbral <- 0.5
clases_predichas <- ifelse(predicciones_prob > umbral, 1, 0)

# Asegurar que ambos sean factores con los mismos niveles para la matriz
test_set$target <- factor(test_set$target, levels = c(0, 1))
clases_predichas <- factor(clases_predichas, levels = c(0, 1))

# 2. Generar la matriz de confusión
matriz_confusion <- table(Real = test_set$target, Predicho = clases_predichas)
print(matriz_confusion)

# 3. Extraer valores
VN <- matriz_confusion[1,1] 
FP <- matriz_confusion[1,2] 
FN <- matriz_confusion[2,1] 
VP <- matriz_confusion[2,2] 

# 4. Calcular métricas
eficiencia  <- VP / (VP + FN)
pureza        <- VP / (VP + FP)
especificidad <- VN / (VN + FP)
producto <- pureza * eficiencia
accuracy  <- (VP + VN) / (VP + VN + FP + FN)

cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n") 


#CURVA ROC
roc_obj <- roc(test_set$target, predicciones_prob)


roc_data <- data.frame(
  TPR = roc_obj$sensitivities,
  FPR = 1 - roc_obj$specificities
)


grafico<-ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#008B8B", size = 1.2) +
  geom_area(fill = "#AFEEEE", alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.65, y = 0.35, 
           label = paste("AUC == ", round(auc(roc_obj), 4)), 
           parse = TRUE, size = 6, color = "#008B8B", fontface = "bold") +
  labs(title = "ROC Curve - XGBoost",
       x = "False Positive Rate (1-Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("roc_xgb.png", plot = grafico, width = 8, height = 6, dpi = 300)


cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n") 

cat("El AUC es:", round(auc(roc_obj), 4), "\n")

