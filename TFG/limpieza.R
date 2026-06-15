# ============================================================
#  LIMPIEZA Y ANÁLISIS EXPLORATORIO
# ============================================================

datos <- read.csv("datos.csv", stringsAsFactors = FALSE)
cat("Filas iniciales:", nrow(datos), "\n\n")

variables_num <- names(datos)[sapply(datos, is.numeric)]
variables_ana <- c("slc_nu_score",
                   "slc_pfp_track_score",
                   "slc_pfp_track_chi2_muon",
                   "slc_pfp_track_chi2_proton",
                   "slc_pfp_shower_length",
                   "slc_pfp_shower_opening_angle",
                   "slc_pfp_shower_energy",
                   "slc_pfp_shower_dEdx")

# ============================================================
# 1. HISTOGRAMAS CON DATOS SIN LIMPIAR
# ============================================================

par(mar = c(4, 4, 3, 1))

for (v in variables_ana) {
  vals <- datos[[v]]
  vals <- vals[!is.na(vals) & is.finite(vals)]
  hist(vals,
       breaks = 80,
       col    = "#AFEEEE",
       border = "white",
       main   = v,
       xlab   = v,
       ylab   = "Frecuencia")
  abline(v = mean(vals),   col = "red",    lwd = 2, lty = 2)
  abline(v = median(vals), col = "orange", lwd = 2, lty = 2)
}
# ============================================================
# 2. MEDIA Y DESVIACIÓN ANTES DE LIMPIAR
# ============================================================
cat("--- Estadísticos ANTES de limpiar ---\n")

resumen_antes <- data.frame(
  variable = variables_ana,
  media    = sapply(variables_ana, function(v) {
    vals <- datos[[v]][is.finite(datos[[v]])]
    format(mean(vals, na.rm = TRUE), scientific = FALSE, digits = 4)
  }),
  desv = sapply(variables_ana, function(v) {
    vals <- datos[[v]][is.finite(datos[[v]])]
    format(sd(vals, na.rm = TRUE), scientific = FALSE, digits = 4)
  })
)
print(resumen_antes, row.names = FALSE)
cat("\n")

#---------------------------------------------


hist(datos$slc_pfp_shower_opening_angle,
     breaks = 10,
     col    = "#AFEEEE",
     border = "white",
     main   = "slc_pfp_shower_opening_angle",
     xlab   = "Ángulo de apertura (rad)",
     ylab   = "Frecuencia",
     ylim   = c(0, 1000),
     xlim   = c(-1050, 100))

legend("topright",
       legend = sprintf("-999 (n = %d)", sum(datos$slc_pfp_shower_opening_angle == -999, na.rm = TRUE)),
       bty = "n")
# ============================================================
# 3. LIMPIEZA: VALORES CENTINELA -999
# ============================================================

cat("-999 por variable:\n")
for (v in variables_num) {
  n <- sum(datos[[v]] == -999, na.rm = TRUE)
  if (n > 0) {
    cat(sprintf("  %-40s %d (%.2f%%)\n", v, n, 100 * n / nrow(datos)))
  }
}
cat("\n")

for (v in variables_num) {
  datos[[v]][datos[[v]] == -999] <- NA
}

filas_antes  <- nrow(datos)
datos_limpios <- datos[complete.cases(datos[, variables_num]), ]
filas_elim   <- filas_antes - nrow(datos_limpios)

cat(sprintf("Filas eliminadas por -999: %d (%.2f%%)\n",
            filas_elim, 100 * filas_elim / filas_antes))
cat("Filas tras limpieza centinelas:", nrow(datos_limpios), "\n\n")

# ============================================================
# 4. MEDIA Y DESVIACIÓN TRAS LIMPIEZA
# ============================================================

cat("--- Estadísticos tras limpieza de centinelas ---\n")
resumen_centinelas <- data.frame(
  variable = variables_ana,
  media    = sapply(variables_ana, function(v)
    format(mean(datos_limpios[[v]], na.rm = TRUE), scientific = FALSE, digits = 4)),
  desv     = sapply(variables_ana, function(v)
    format(sd(datos_limpios[[v]],   na.rm = TRUE), scientific = FALSE, digits = 4))
)
print(resumen_centinelas, row.names = FALSE)
cat("\n")




#HISTOGRAMA SHOWER DEDX


hist(datos_limpios$slc_pfp_shower_dEdx,
     breaks = 500,
     col    = "#AFEEEE",
     border = "white",
     main   = "slc_pfp_shower_dEdx",
     xlab   = "dE/dx (MeV/cm)",
     ylab   = "Frecuencia",
     ylim=c(0,100),
     xlim   = c(0, 400000))
abline(v = mean(datos_limpios$slc_pfp_shower_dEdx,   na.rm = TRUE), col = "red",    lwd = 2, lty = 2)
abline(v = median(datos_limpios$slc_pfp_shower_dEdx, na.rm = TRUE), col = "orange", lwd = 2, lty = 2)
legend("topright",
       legend = c(sprintf("Media   = %.1f", mean(datos_limpios$slc_pfp_shower_dEdx,   na.rm = TRUE)),
                  sprintf("Mediana = %.2f", median(datos_limpios$slc_pfp_shower_dEdx, na.rm = TRUE))),
       col = c("red", "orange"), lwd = 2, lty = 2, bty = "n")


hist(datos_limpios$slc_pfp_shower_dEdx,
     breaks = 200,
     col    = "#AFEEEE",
     border = "white",
     main   = "slc_pfp_shower_dEdx",
     xlab   = "dE/dx (MeV/cm)",
     ylab   = "Frecuencia")
abline(v = 20, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = sprintf("Corte dEdx = 20  (outliers n = %d)", sum(datos_limpios$slc_pfp_shower_dEdx > 20)),
       col = "red", lwd = 2, lty = 2, bty = "n")



# 5. LIMPIEZA: OUTLIERS dEdx 
# ============================================================

# Corte dEdx
n_antes_dedx <- nrow(datos_limpios)
datos_limpios <- datos_limpios[datos_limpios$slc_pfp_shower_dEdx < 20, ]

cat(sprintf("Filas eliminadas por dEdx >= 20: %d (%.2f%%)\n",
            n_antes_dedx - nrow(datos_limpios),
            100 * (n_antes_dedx - nrow(datos_limpios)) / n_antes_dedx))

# ============================================================
# 6. MEDIA Y DESVIACIÓN TRAS LIMPIEZA DE OUTLIERS
# ============================================================

cat("--- Estadísticos DESPUÉS de limpiar outliers ---\n")
resumen_despues <- data.frame(
  variable = variables_ana,
  media    = sapply(variables_ana, function(v)
    format(mean(datos_limpios[[v]], na.rm = TRUE), scientific = FALSE, digits = 4)),
  desv     = sapply(variables_ana, function(v)
    format(sd(datos_limpios[[v]],   na.rm = TRUE), scientific = FALSE, digits = 4))
)
print(resumen_despues, row.names = FALSE)
cat("\n")

# ============================================================
# 7. GUARDAR DATASET LIMPIO
# ============================================================

write.csv(datos_limpios, "datos_limpios.csv", row.names = FALSE)
