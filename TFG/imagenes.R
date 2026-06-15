library(ggplot2)
library(MASS)   

set.seed(123)

#-------------LDA---------------------


# --- Datos ---
n1 <- 200; n2 <- 200

mu1 <- c(2.2, 3.0)
mu2 <- c(5.0, 6.3)
Sigma <- matrix(c(0.6, 0.15, 0.15, 0.5), 2, 2) #misma matriz de covarianza

#distribucion normal multivariante
X1 <- mvrnorm(n1, mu1, Sigma)
X2 <- mvrnorm(n2, mu2, Sigma)

df <- data.frame(
  x1    = c(X1[,1], X2[,1]),
  x2    = c(X1[,2], X2[,2]),
  grupo = factor(rep(c("Clase 1", "Clase 2"), c(n1, n2)))
)

# --- Fisher LDA ---
m1 <- colMeans(X1)
m2 <- colMeans(X2)
W  <- cov(X1) * (n1 - 1) + cov(X2) * (n2 - 1) #matriz intra-clase 

a <- solve(W) %*% (m2 - m1) #dir maxima separacion
a <- a / sqrt(sum(a^2))

M <- (m1 + m2) / 2 #punto medio entre clases
b <- c(-a[2], a[1]) #perpendicular a a

t_vals      <- seq(-9, 9, length.out = 300)
proj_df     <- data.frame(x = M[1] + t_vals * a[1], y = M[2] + t_vals * a[2])
boundary_df <- data.frame(x = M[1] + t_vals * b[1], y = M[2] + t_vals * b[2])

# --- Puntos especiales para leyenda ---
special_df <- data.frame(
  x     = c(m1[1],     m2[1],     M[1]),
  y     = c(m1[2],     m2[2],     M[2]),
  grupo = c("Media μ_1", "Media μ_2", "Punto Medio (M)")
)

fill_vals <- c(
  "Clase 1"          = "#20B2AA",
  "Clase 2"          = "#EE44EE",
  "Media μ_1"         = "#008B8B",
  "Media μ_2"         = "#CC00CC",
  "Punto Medio (M)"  = "black"
)

# --- Gráfico ---
ggplot() +
  geom_line(data = proj_df,
            aes(x = x, y = y, linetype = "Eje de Proyección (a)"),
            color = "gray40", linewidth = 0.9) +
  geom_line(data = boundary_df,
            aes(x = x, y = y, color = "Frontera de Decisión (y(x) = 0)"),
            linewidth = 1.3) +
  geom_point(data = df,
             aes(x = x1, y = x2, fill = grupo),
             shape = 21, size = 3, color = "white", alpha = 0.85) +
  geom_point(data = special_df,
             aes(x = x, y = y, fill = grupo),
             shape = 21, size = 4.5, color = "white") +

  scale_fill_manual(name = NULL, values = fill_vals) +
  scale_color_manual(
    name   = NULL,
    values = c("Frontera de Decisión (y(x) = 0)" = "lightgreen")
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Eje de Proyección (a)" = "dashed")
  ) +

  coord_cartesian(xlim = c(-2, 10), ylim = c(0, 10)) +
  scale_x_continuous(breaks = seq(-2, 10, 2)) +
  scale_y_continuous(breaks = seq(0,  10, 2)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major  = element_line(color = "gray88"),
    panel.grid.minor  = element_blank(),
    legend.position   = c(0.19, 0.80),
    legend.background = element_rect(fill = "white", color = "gray70"),
    legend.key        = element_rect(fill = "white"),
    legend.text       = element_text(size = 7.5),
    legend.key.size   = unit(3.5, "mm"),
    legend.spacing.y  = unit(0.5, "mm"),
    plot.title        = element_blank()
  ) +
  labs(x = expression(X[1]), y = expression(X[2]))




#-----------QDA-----------------


n <- 300

# Clase cyan
x_c <- c(rnorm(200, -0.8, 0.7), rnorm(100, 1.3, 0.9))
y_c <- c(rnorm(200,  1.1, 0.6), rnorm(100, -0.2, 0.6))

# Clase magenta
x_m <- rnorm(n, -1.1, 0.8)
y_m <- rnorm(n, -1.0, 0.6)

df <- data.frame(
  x     = c(x_c, x_m),
  y     = c(y_c, y_m),
  clase = factor(rep(c("cyan", "magenta"), c(length(x_c), length(x_m))))
)

# Modelo QDA 
modelo_qda <- qda(clase ~ x + y, data = df)


grid <- expand.grid(
  x = seq(-3, 3.5, by = 0.04),
  y = seq(-3.2, 3.6, by = 0.04)
)

# Probabilidad posterior de clase magenta
pred       <- predict(modelo_qda, newdata = grid)
grid$prob  <- pred$posterior[, "magenta"]


ggplot() +
  geom_raster(data = grid,
              aes(x = x, y = y, fill = prob),
              interpolate = TRUE) +
  scale_fill_gradientn(
    colors = c("#40E0D0", "white", "#EE82EE"),
    values = c(0, 0.5, 1)
  ) +
  geom_point(data = df,
             aes(x = x, y = y, color = clase),
             shape = 1, size = 2.5, stroke = 1.2) +
  scale_color_manual(values = c("cyan" = "#20B2AA", "magenta" = "#EE44EE")) +

  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-3, 3.5), ylim = c(-3.2, 3.6)) +
  labs(x = NULL, y = NULL)






#------------FUNCION SIGMOIDE LOGIT-------------------


# Función sigmoide
sigmoid <- function(x, b0, b1) 1 / (1 + exp(-(b0 + b1 * x)))

x <- seq(-10, 10, length.out = 500)


df <- rbind(
  data.frame(x = x, y = sigmoid(x, 0,  1), curva = "c1"),
  data.frame(x = x, y = sigmoid(x, 0,  3), curva = "c2"),
  data.frame(x = x, y = sigmoid(x, -4, 1), curva = "c3")
)

df$curva <- factor(df$curva, levels = c("c1", "c2", "c3"))


ggplot(df, aes(x = x, y = y, color = curva, linetype = curva)) +
  geom_hline(yintercept = c(0, 0.5, 1), color = "gray70", linewidth = 0.4, linetype = "dotted") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    name   = NULL,
    values = c("c1" = "blue", "c2" = "red", "c3" = "darkgreen"),
    labels = expression(
      alpha==0 ~ "," ~ beta==1,
      alpha==0 ~ "," ~ beta==3 ~ "(más rápida)",
      alpha==-4 ~ "," ~ beta==1 ~ "(desplazada)"
    )
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("c1" = "solid", "c2" = "dashed", "c3" = "dotted"),
    labels = expression(
      alpha==0 ~ "," ~ beta==1,
      alpha==0 ~ "," ~ beta==3 ~ "(más rápida)",
      alpha==-4 ~ "," ~ beta==1 ~ "(desplazada)"
    )
  ) +
  scale_x_continuous(breaks = seq(-10, 10, 2.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major  = element_line(color = "gray88"),
    panel.grid.minor  = element_blank(),
    legend.position   = c(0.22, 0.83),
    legend.background = element_rect(fill = "white", color = "gray70"),
    legend.key        = element_rect(fill = "white"),
    legend.text       = element_text(size = 11),
    legend.key.size   = unit(8, "mm"),
    plot.title        = element_text(hjust = 0.5, size = 11)
  ) +
  labs(
    title = "Forma de la Probabilidad en el Modelo Logit (Función Sigmoide)",
    x     = "Variable Explicativa (X)",
    y     = "Probabilidad P(Y=1)"
  )
