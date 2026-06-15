## TEST DE SEPARACIÓN 

library(tidyverse)

# Buscamos variables prometedoras con ayuda de tests

datos <- read_csv("datos.csv") %>%
  mutate(track_score_std = scale(slc_pfp_track_score)) #normalizamos

#--- HISTOGRAMA------

ggplot(datos, aes(x = track_score_std, fill = factor(target))) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("0" = "#40E0D0", "1" = "#FF6F61"),
                    labels = c("0" = "Fondo", "1" = "ν_e CC")) +
  labs(title = "Track score estandarizado",
       x = "Track score", y = "Frecuencia", fill = "Tipo") +
  theme_minimal()

# --- PODER DE SEPARACIÓN (η) ---
stats <- datos %>% 
  group_by(target) %>% #agrupa los datos por clase 
  summarise(mu = mean(track_score_std, na.rm = TRUE),
            v  = var(track_score_std,  na.rm = TRUE))

# Aplicar la fórmula del Poder de Separación
# η = |media1 - media2| / raiz(varianza1 + varianza2)

eta <- with(stats, abs(diff(mu)) / sqrt(sum(v)))
print(paste("Poder de separación η =", round(eta, 4)))


#TODAS LAS VARIABLES

variables_reco <- datos %>%
  select(starts_with("slc_")) %>%
  names()

resultados <- map_dfr(variables_reco, function(var) {
  stats <- datos %>%
    group_by(target) %>%
    summarise(mu = mean(.data[[var]], na.rm = TRUE),
              v  = var(.data[[var]],  na.rm = TRUE))
  
  tibble(
    variable = var,
    eta = with(stats, abs(diff(mu)) / sqrt(sum(v)))
  )
}) %>%
  arrange(desc(eta))

print(resultados)



#TEST CHI^2 DE PEARSON 

signal <- datos$slc_pfp_track_score[datos$target == 1]
fondo  <- datos$slc_pfp_track_score[datos$target == 0]

breaks <- seq(0, 1, length.out = 51)

n_signal <- hist(signal, breaks = breaks, plot = FALSE)$counts
n_fondo  <- hist(fondo,  breaks = breaks, plot = FALSE)$counts

# Eliminar bins donde ambos son 0
mask <- (n_signal + n_fondo) > 0
n_signal <- n_signal[mask]
n_fondo  <- n_fondo[mask]

chisq.test(rbind(n_signal, n_fondo))





signal <- datos$slc_pfp_shower_dEdx[datos$target == 1]
fondo  <- datos$slc_pfp_shower_dEdx[datos$target == 0]

breaks <- seq(-1, 20, length.out = 51)

# Filtrar también los valores fuera del rango
signal <- signal[signal >= -1 & signal <= 20]
fondo  <- fondo[fondo >= -1 & fondo <= 20]

n_signal <- hist(signal, breaks = breaks, plot = FALSE)$counts
n_fondo  <- hist(fondo,  breaks = breaks, plot = FALSE)$counts

mask <- (n_signal + n_fondo) > 0
n_signal <- n_signal[mask]
n_fondo  <- n_fondo[mask]

chisq.test(rbind(n_signal, n_fondo), simulate.p.value = TRUE, B = 10000)

#TODAS LAS VARIABLES


variables_slc <- names(datos)[grepl("^slc", names(datos))]

resultados_chi2 <- data.frame(variable = variables_slc, chi2 = NA, p_valor = NA)

for (i in seq_along(variables_slc)) {
  
  signal <- datos[[variables_slc[i]]][datos$target == 1]
  fondo  <- datos[[variables_slc[i]]][datos$target == 0]
  
  signal <- signal[is.finite(signal)]
  fondo  <- fondo[is.finite(fondo)]
  
  breaks <- seq(min(datos[[variables_slc[i]]], na.rm = TRUE) - 0.01,
                max(datos[[variables_slc[i]]], na.rm = TRUE) + 0.01,
                length.out = 51)
  
  n_signal <- hist(signal, breaks = breaks, plot = FALSE)$counts
  n_fondo  <- hist(fondo,  breaks = breaks, plot = FALSE)$counts
  
  mask <- (n_signal + n_fondo) > 0
  n_signal <- n_signal[mask]
  n_fondo  <- n_fondo[mask]
  
  test <- chisq.test(rbind(n_signal, n_fondo), simulate.p.value = TRUE, B = 10000)
  
  resultados_chi2$chi2[i]    <- round(test$statistic, 2)
  resultados_chi2$p_valor[i] <- test$p.value
}

resultados_chi2[order(resultados_chi2$p_valor), ]


#K-S

signal <- datos$slc_pfp_shower_dEdx[datos$target == 1]
fondo  <- datos$slc_pfp_shower_dEdx[datos$target == 0]

# Eliminar NA e Inf
signal <- signal[is.finite(signal)]
fondo  <- fondo[is.finite(fondo)]

ks.test(signal, fondo)

#TODAS VARIABLES

variables_slc <- names(datos)[grepl("^slc", names(datos))]

resultados_ks <- data.frame(variable = variables_slc, D = NA, p_valor = NA)

for (i in seq_along(variables_slc)) {
  
  signal <- datos[[variables_slc[i]]][datos$target == 1]
  fondo  <- datos[[variables_slc[i]]][datos$target == 0]
  
  signal <- signal[is.finite(signal)]
  fondo  <- fondo[is.finite(fondo)]
  
  test <- ks.test(signal, fondo)
  
  resultados_ks$D[i]       <- round(test$statistic, 4)
  resultados_ks$p_valor[i] <- test$p.value
}
resultados_ks$p_valor <- ifelse(
  resultados_ks$p_valor == 0,
  "< 2.2e-16",
  formatC(resultados_ks$p_valor, format = "e", digits = 4)
)

resultados_ks[order(resultados_ks$D, decreasing = TRUE), ]




