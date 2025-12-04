# Librerías para gráficos
library(ggplot2)
library(latex2exp)
library(patchwork)

set.seed(72) # Libreta de Delfina Kiss, 72/23

# ----- Ejercicio 1.5: simulación Monte Carlo -----
theta_real <- 0.25
alpha <- 0.05
N <- 5000

n_values <- c(3, 10, 30, 50, 100)

coberturas <- numeric(length(n_values))

for (i in 1:length(n_values)) {
  n <- n_values[i]
  t_por_sim <- rbinom(n = N, size = n, prob = theta_real)
  theta_hat <- t_por_sim / n
  se <- sqrt((theta_hat * (1 - theta_hat)) / n)
  z_alpha <- qnorm(1 - alpha / 2)
  lim_inf <- theta_hat - z_alpha * se
  lim_sup <- theta_hat + z_alpha * se
  
  aciertos <- (theta_real >= lim_inf) & (theta_real <= lim_sup)
  coberturas[i] <- mean(aciertos)
}

resultados <- data.frame(n = n_values, Cobertura = coberturas)
print(resultados)

rm(list = ls()) # Borramos variables de ambiente para el próximo ejercicio

# ----- Ejercicio 2.3: gráficos de p -----

se_real <- 0.9
sp_real <- 0.95
theta_real <- 0.25

# Definimos la función p
p_func <- function(se, sp, theta) {
  return(theta * (se + sp - 1) + 1 - sp)
}

p_actual <- p_func(se_real, sp_real, theta_real)

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
y_se <- p_func(x_se, sp_real, theta_real)
p1 <- crear_grafico_individual(x_se, y_se, "$Se$", "Variación en $Se$", se_real)

# Segundo gráfico variando Sp
x_sp <- seq(0, 1, length.out = 100)
y_sp <- p_func(se_real, x_sp, theta_real)
p2 <- crear_grafico_individual(x_sp, y_sp, "$Sp$", "Variación en $Sp$", sp_real, show_yaxis = FALSE)

# Tercer gráfico variando theta
x_th <- seq(0, 1, length.out = 100)
y_th <- p_func(se_real, sp_real, x_th)
p3 <- crear_grafico_individual(x_th, y_th, "$\\theta$", "Variación en $\\theta$", theta_real, show_yaxis = FALSE)

# Combinamos los tres gráficos
grafico_final <- p1 + p2 + p3 + 
  plot_annotation(
    title = TeX("Análisis de Sensibilidad de $p$"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, family = "serif")
      )
  )

print(grafico_final)



