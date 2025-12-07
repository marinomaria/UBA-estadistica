# Librerías para gráficos
library(ggplot2)
library(latex2exp)
library(patchwork)

set.seed(72) # Libreta de Delfina Kiss, 72/23
# Fijamos los parámetros de Theta, Se y Sp del enunciado
THETA0 <- 0.25
SE0 <- 0.90
SP0 <- 0.95

# ----- Ejercicio 1.5: simulación Monte Carlo -----
alpha <- 0.05
N <- 5000

n_values <- c(3, 10, 30, 50, 100)

coberturas <- numeric(length(n_values))

for (i in 1:length(n_values)) {
  n <- n_values[i]
  t_por_sim <- rbinom(n = N, size = n, prob = THETA0)
  theta_hat <- t_por_sim / n
  se <- sqrt((theta_hat * (1 - theta_hat)) / n)
  z_alpha <- qnorm(1 - alpha / 2)
  lim_inf <- theta_hat - z_alpha * se
  lim_sup <- theta_hat + z_alpha * se
  
  aciertos <- (THETA0 >= lim_inf) & (THETA0 <= lim_sup)
  coberturas[i] <- mean(aciertos)
}

resultados <- data.frame(n = n_values, Cobertura = coberturas)
print(resultados)

# ----- Ejercicio 2.0.3: gráficos de p -----

# Definimos la función p
p_func <- function(se, sp, theta) {
  return(theta * (se + sp - 1) + 1 - sp)
}

p_actual <- p_func(SE0, SP0, THETA0)

crear_grafico_individual <- function(
    x_range,
    y_vals,
    x_lab, 
    title,
    intercept_x,
    show_yaxis = TRUE
    ) {
  
  data <- data.frame(x = x_range, y = y_vals)
  
  plot_obj <- ggplot(data, aes(x = x, y = y)) +
    geom_line(color = "#78911d", size = 1) + 
    # Agregamos un punto con el valor real del param
    geom_point(aes(x = intercept_x, y = p_actual), 
               color = "#E74C3C", size = 1.5) +
    labs(
      title = TeX(title),
      x = TeX(x_lab),
      y = if(show_yaxis) TeX("$p(Se, Sp, \\theta)$") else NULL
    ) +
    theme_bw(base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11)
    ) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE)
  
  if (!show_yaxis) {
    plot_obj <- plot_obj + 
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  return(plot_obj)
}

# Primer gráfico variando Se
x_se <- seq(0, 1, length.out = 100)
y_se <- p_func(x_se, SP0, THETA0)
p1 <- crear_grafico_individual(x_se, y_se, "$Se$", "Variación en $Se$", SE0)

# Segundo gráfico variando Sp
x_sp <- seq(0, 1, length.out = 100)
y_sp <- p_func(SE0, x_sp, THETA0)
p2 <- crear_grafico_individual(x_sp, y_sp, "$Sp$", "Variación en $Sp$", SP0, show_yaxis = FALSE)

# Tercer gráfico variando theta
x_th <- seq(0, 1, length.out = 100)
y_th <- p_func(SE0, SP0, x_th)
p3 <- crear_grafico_individual(x_th, y_th, "$\\theta$", "Variación en $\\theta$", THETA0, show_yaxis = FALSE)

# Combinamos los tres gráficos
grafico_final <- p1 + p2 + p3 +
  plot_annotation(
    title = TeX("Análisis de sensibilidad de $p(Se, Sp, \\theta)$"),
    theme = theme(
      plot.title = element_text(family = "serif", face = "bold", size = 16, hjust = 0.5)
    )
  )

print(grafico_final)

# ----- Ejercicio 2.1.6: comparación de ECM -----

ecm_emv_perfecto <- function(n, theta) {
  return(theta * (1 - theta) / n)
}

ecm_mom_imperfecto <- function(n, theta, se, sp) {
  p <- theta * (se + sp - 1) + 1 - sp
  return((1 / n) * (p * (1 - p) / (se + sp - 1)^2))
}

rango_n <- seq(1, 100, by = 1)

matriz_de_comparacion <- data.frame(
  n = rep(rango_n, 2),
  value = c(
    ecm_emv_perfecto(rango_n, THETA0), 
    ecm_mom_imperfecto(rango_n, THETA0, SE0, SP0)
  ),
  func_type = factor(rep(c("ECM EMV", "ECM MoM"), each = length(rango_n)))
)

ggplot(matriz_de_comparacion, aes(x = n, y = value, color = func_type)) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("ECM EMV" = "#78911d", "ECM MoM" = "#E74C3C"),
    labels = c("ECM EMV" = TeX("$ECM_{EMV}$ - Test Perfecto"), 
               "ECM MoM" = TeX("$ECM_{MoM}$ - Test Imperfecto")),
    name = NULL
  ) +

  labs(title = "Comparación de convergencia de ECM",
    x = TeX("$n$"),
    y = TeX("$ECM(n)$")
  ) +
  
  theme_bw(base_family = "serif") +
  theme(
    legend.position = c(0.75, 0.8),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    panel.grid.minor = element_blank()
  ) +
  
  coord_cartesian(expand = FALSE, ylim = c(0, 0.3), xlim = c(1, 50))

# ----- Ejercicio 2.1.7: simulación Monte Carlo -----

n <- 100
Nrep <- 1000

# Usamos la fórmula derivada en el inciso 2.0.2
p_val <- THETA0 * (SE0 + SP0 - 1) + (1 - SP0)

T_sum <- rbinom(n = Nrep, size = n, prob = p_val)

T_promedio <- T_sum / n

theta_mom <- (T_promedio + SP0 - 1) / (SE0 + SP0 - 1)

sesgo_teorico <- 0
sesgo_empirico <- mean(theta_mom) - THETA0
cat("Sesgo Teórico:", sesgo_teorico, "\t Sesgo Empírico:", sesgo_empirico, "\n")

var_teorica <- (p_val * (1 - p_val)) / (n * (SE0 + SP0 - 1)^2)
var_empirica <- var(theta_mom)
cat("Varianza Teórica:", var_teorica, "\t Varianza Empírica:", var_empirica, "\n")

# ----- Ejercicio 2.1.8: bootstrap -----

n <- 10 # Por enunciado
B <- 1000

# Calculamos p teórico
p_val <- p_func(SE0, SP0, THETA0)

# Generamos una muestra aleatoria Bi(n, p teórico)
muestra_base <- rbinom(n = n, size = 1, prob = p_val)

theta_mom_b <- numeric(B)

for(i in 1:B) {
  # Remuestreo con reposición de la muestra ORIGINAL
  muestra_b <- sample(muestra_base, size = n, replace = TRUE)
  
  # Proporción de positivos en la muestra bootstrap
  p_hat_b <- mean(muestra_b)
  
  theta_mom_b[i] <- (p_hat_b + SP0 - 1) / (SE0 + SP0 - 1)
}

datos_bootstrap <- data.frame(theta_est = theta_mom_b)

# Graficamos para ver la distribución del estimador
ggplot(datos_bootstrap, aes(x = theta_est)) +
  geom_histogram(bins = 15, fill = "#78911d",  alpha = 0.7) +
  geom_vline(xintercept = THETA0, color = "#E74C3C", linetype = "dashed", size = 1) +
  labs(title = "Distribución Bootstrap del Estimador de Momentos",
       x = TeX("$\\hat{\\theta}_{MoM}$"),
       y = "Frecuencia") +
  theme_bw(base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11))


# ----- Ejercicio 3.1.2: obtención de datos para evaluar test -----

n <- 100
theta_pre <- 0.2
theta_post <- 0.15
p_pre <- p_func(SE0, SP0, theta_pre)
p_post <- p_func(SE0, SP0, theta_post)

X_pre <- rbinom(n = n, size = 1, p_pre)
X_post <- rbinom(n = n, size = 1, p_post)

prom_pre <- mean(X_pre)
prom_post <- mean(X_post)

# Bajo H0 asumimos que p_pre = p_post, así que juntamos todos los positivos
p_hat_H0 <- (sum(X_pre) + sum(X_post)) / (2 * n)

res <- data.frame(
  Pre  = c(p_pre, prom_pre), 
  Post = c(p_post, prom_post)
)
row.names(res) <- c("p real", "Promedio Muestral")
print(res)
cat(sprintf("Estimador combinado de p (H0): %.4f\n", p_hat_H0))

# ----- Ejercicio 3.1.3: nivel empírico del test -----

n_scenarios <- matrix(c(
  15, 20,
  56, 40,
  100, 100, 
  10, 100, 
  500, 458,
  1245, 1455
), ncol = 2, byrow = TRUE)

prop_rechazos <- numeric(nrow(n_scenarios))
alpha <- 0.05
z_alpha <- qnorm(1 - alpha/2) 

for (i in 1:nrow(n_scenarios)) {
  n_pre <- n_scenarios[i, 1]
  n_post <- n_scenarios[i, 2]
  
  x_pre <- rbinom(n = Nrep, size = n_pre, prob = p_func(SE0, SP0, theta_pre))
  x_post <- rbinom(n = Nrep, size = n_post, prob = p_func(SE0, SP0, theta_post))
  
  p_hat_pre <- x_pre / n_pre
  p_hat_post <- x_post / n_post
  
  p_pool <- (x_pre + x_post) / (n_pre + n_post)
  
  pivote_h0 <- (p_hat_post - p_hat_pre) / sqrt(p_pool * (1 - p_pool) * (1/n_pre + 1/n_post))
  
  rechazos <- abs(pivote_h0) >= z_alpha
  
  prop_rechazos[i] <- mean(rechazos)
}

print(data.frame(n_pre = n_scenarios[,1], n_post = n_scenarios[,2], Potencia = prop_rechazos))






