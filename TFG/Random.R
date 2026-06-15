# RANDOM FOREST 


library(caTools)
library(randomForest)
library(dplyr)
library(readr)
library(pROC)


datos <- read.csv("datos_limpios.csv")

datos$target <- as.factor(datos$target)


set.seed(123) 

division <- sample.split(datos$target, SplitRatio = 0.70)

train_set <- subset(datos, division == TRUE)
test_set <- subset(datos, division == FALSE)


modelo_rf <- randomForest(target ~ slc_pfp_track_score +
                            slc_pfp_shower_energy +
                            slc_pfp_track_chi2_muon +
                            slc_nu_score+
                            slc_pfp_shower_length +
                            slc_pfp_shower_opening_angle,
                          data = train_set, 
                          ntree = 100)

print(modelo_rf)

predicciones <- predict(modelo_rf, test_set)

matriz<- table(Real = test_set$target, Predicho = predicciones)
print(matriz)

VN <- matriz[1,1] 
FP <- matriz[1,2] 
FN <- matriz[2,1] 
VP <- matriz[2,2] 


eficiencia <- VP / (VP + FN)
especificidad <- VN / (VN + FP)
pureza <- VP / (VP + FP)
producto<- eficiencia * pureza
accuracy <- (VP + VN) / (VP + VN + FP + FN)



probabilidades <- predict(modelo_rf, test_set, type = "prob")[,2]


#CURVA ROC
roc_obj <- roc(test_set$target, probabilidades)


roc_data <- data.frame(
  TPR = roc_obj$sensitivities,
  FPR = 1 - roc_obj$specificities
)

grafico<-ggplot(roc_data, aes(x = FPR, y = TPR)) +
  
  geom_line(color = "#008B8B", size = 1) +

  geom_area(fill = "#AFEEEE", alpha = 0.4) +
  # Línea de azar
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +

  annotate("text", x = 0.65, y = 0.35, 
           label = paste("AUC == ", round(auc(roc_obj), 4)), 
           parse = TRUE, size = 6, color = "#008B8B", fontface = "bold") +
  labs(title = "ROC Curve - Random Forest",
       x = "False Positive Rate (1- Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray92")
   
  )
ggsave("roc_rf.png", plot = grafico, width = 8, height = 6, dpi = 300)


cat("Eficiencia:", eficiencia, "\n")
cat("Especificidad:", especificidad, "\n")
cat("Pureza:", pureza, "\n")
cat("Ef * Pureza:", producto, "\n")
cat("Accuracy:", accuracy, "\n") 


cat("El AUC es:", round(auc(roc_obj), 4), "\n")


