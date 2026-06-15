##ANALISIS EXPLORATORIO 

library(ggplot2)

datos<-read.csv("datos.csv")

#HISTOGRAMAS

ggplot(datos, aes(x = slc_pfp_track_score, fill = factor(target))) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("0" = "#40E0D0", "1" = "#FF6F61"),
                    labels = c("0" = "Fondo", "1" = "Señal")) +
  labs(x = "Track Score", y = "Frecuencia", fill = "") +
  theme_minimal()

ggsave("track_score.png", 
       width = 8, height = 5, units = "in", dpi = 300)


ggplot(datos, aes(x = slc_pfp_shower_energy, fill = factor(target))) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("0" = "#40E0D0", "1" = "#FF6F61"),
                    labels = c("0" = "Fondo", "1" = "Señal")) +
  scale_x_continuous(limits = c(0, 3000)) +
  labs(x = "Shower Energy", y = "Frecuencia", fill = "") +
  theme_minimal()

ggsave("shower_energy.png", 
       width = 8, height = 5, units = "in", dpi = 300)


ggplot(datos, aes(x = slc_pfp_shower_length, fill = factor(target))) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("0" = "#40E0D0", "1" = "#FF6F61"),
                    labels = c("0" = "Fondo", "1" = "Señal")) +
  scale_x_continuous(limits = c(0, 300)) +
  labs(x = "Shower Length", y = "Frecuencia", fill = "") +
  theme_minimal()


ggsave("shower_length.png", 
       width = 8, height = 5, units = "in", dpi = 300)

#DIAGRAMAS DE DISPERSION

ggplot(datos, aes(x = slc_pfp_shower_opening_angle, y = slc_pfp_shower_length, color = factor(target))) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("0" = "#40E0D0", "1" = "#FF6F61"), 
                     labels = c("Fondo", "Señal")) +
  xlim(-0.1, 1) +
  ylim(0, 300) +
  labs(title = "Scatter plot: opening_angle vs shower_length",
       x = "Shower Opening Angle", y = "Shower Length", color = "Tipo de suceso") +
  theme_minimal()

ggsave("correlacion.png", 
       width = 8, height = 5, units = "in", dpi = 300)


library(GGally)

variables_reco <- c("slc_pfp_track_score", "slc_pfp_shower_energy",
                    "slc_nu_score", "slc_pfp_track_chi2_muon",
                    "slc_pfp_shower_length", "slc_pfp_track_chi2_proton",
                    "slc_pfp_shower_dEdx")

variables_corr <- c("slc_pfp_track_score", "slc_pfp_shower_energy",
                    "slc_nu_score", "slc_pfp_track_chi2_muon",
                    "slc_pfp_shower_length", "slc_pfp_track_chi2_proton",
                    "slc_pfp_shower_dEdx", "slc_pfp_shower_opening_angle")

library(GGally)

ggpairs(datos[, variables_reco],
        aes(color = factor(datos$target), alpha = 0.3)) +
  theme_minimal()

library(corrplot)
corrplot(cor(datos[, variables_corr]),
         method = "color",
         addCoef.col = "black",
         tl.cex = 0.8)



png("matriz_dispersion.png", width = 1200, height = 1200)

colores <- c("0" = "#40E0D0", "1" = "salmon")
colores_por_fila <- colores[as.character(datos$target)]

# Crear datos con límites aplicados
datos_lim <- datos
datos_lim$slc_pfp_track_chi2_proton[datos_lim$slc_pfp_track_chi2_proton > 500] <- NA
datos_lim$slc_pfp_shower_dEdx[datos_lim$slc_pfp_shower_dEdx > 10] <- NA

pairs(
  datos_lim[, variables_reco],
  col = colores_por_fila,
  pch = 19,               
  cex = 0.3,              
  main = "Matriz de Dispersión",
  lower.panel = NULL
)

legend(
  "bottomleft", 
  legend = c("Fondo", "Señal"), 
  pch = 19, 
  col = c("cyan", "salmon"), 
  cex = 0.8,
  bty = "n"
)

dev.off()

summary(datos$slc_pfp_track_chi2_proton)
summary(datos$slc_pfp_shower_dEdx)