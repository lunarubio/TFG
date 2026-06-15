library(ggplot2)
library(tidyverse)
library(dplyr)
library(MASS)
library(pROC)
library(caTools)

datos<-read.csv("datos_limpios.csv")

set.seed(123)

#asegurar que es factor
datos$target <- as.factor(datos$target)


#MODELO LOGIT


# 1. División del dataset (70% entrenamiento, 30% prueba)

division <- sample.split(datos$target, SplitRatio = 0.70)

datos_entrenamiento <- subset(datos, division == TRUE)
datos_prueba        <- subset(datos, division == FALSE)

# Normalizar variables (media 0, SD 1) usando parámetros del entrenamiento
vars_pred <- c("slc_pfp_track_score", "slc_pfp_shower_energy",
               "slc_nu_score", "slc_pfp_track_chi2_muon",
               "slc_pfp_shower_length", "slc_pfp_shower_opening_angle")

medias <- sapply(datos_entrenamiento[, vars_pred], mean)
desv   <- sapply(datos_entrenamiento[, vars_pred], sd)

datos_entrenamiento[, vars_pred] <- scale(datos_entrenamiento[, vars_pred])
datos_prueba[, vars_pred] <- scale(datos_prueba[, vars_pred],
                                   center = medias, scale = desv)

# 2. Entrenar el modelo con el conjunto de entrenamiento

modelo_logit_train <- glm(target ~ slc_pfp_track_score + slc_pfp_shower_energy +
                            slc_nu_score + slc_pfp_track_chi2_muon +
                            slc_pfp_shower_length + slc_pfp_shower_opening_angle,
                          data = datos_entrenamiento,
                          family = binomial(link = "logit"))

summary(modelo_logit_train)

# 3. Predecir probabilidades y clasificarlas con el umbral 0.5 sobre el test

predicciones_test <- predict(modelo_logit_train, newdata = datos_prueba, type = "response")
clases_test       <- ifelse(predicciones_test > 0.5, 1, 0)

# 4. Matriz de confusión 
matriz <- table(Real = datos_prueba$target, Predicho = clases_test)
print(matriz)


VN <- matriz[1,1]
FP <- matriz[1,2]
FN <- matriz[2,1]
VP <- matriz[2,2]

# 5. Métricas

eficiencia  <- VP / (VP + FN)
especificidad <- VN / (VN + FP)
pureza        <- VP / (VP + FP)
producto      <- eficiencia * pureza
accuracy      <- (VP + VN) / (VP + VN + FP + FN)



# 6. Curva ROC basada estrictamente en los datos de prueba
roc_obj <- roc(datos_prueba$target, predicciones_test)



# 2. Preparar los datos para el formato de ggplot
# El objeto roc de pROC contiene las especificidades y eficienciaes
roc_data <- data.frame(
  TPR = roc_obj$sensitivities,
  FPR = 1 - roc_obj$specificities
)

# 3. Generar el gráfico con tu formato
grafico<-ggplot(roc_data, aes(x = FPR, y = TPR)) +
  # Línea turquesa intenso
  geom_line(color = "#008B8B", size = 1) +
  # Área sombreada turquesa claro
  geom_area(fill = "#AFEEEE", alpha = 0.4) +
  # Línea de azar
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # Anotación del AUC dentro del área
  annotate("text", x = 0.65, y = 0.35, 
           label = paste("AUC == ", round(auc(roc_obj), 4)), 
           parse = TRUE, size = 6, color = "#008B8B", fontface = "bold") +
  labs(title = "ROC Curve - Modelo Logit",
       x = "False Positive Rate (1- Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray92")
  )

ggsave("roc_logit.png", plot = grafico, width = 8, height = 6, dpi = 300)

cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n") 
cat("El AUC es:", round(auc(roc_obj), 4), "\n")


exp(coef(modelo_logit_train))

confint(modelo_logit_train)        # IC de los coeficientes β
exp(confint(modelo_logit_train))   # IC del OR

