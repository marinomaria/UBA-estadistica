# Librerías para gráficos
library(ggplot2)
library(latex2exp)
library(patchwork)

set.seed(72) # Libreta de Maria Delfina Kiss, 72/23
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

# ----- Ejercicio 2.1.7: simulación Monte Carlo-----

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



# ----- Ejercicio 2.2: intervalos de confianza bootstrap -----

# Definimos los distintos n
n=c(5,10,20,50,100,1000)

# Método no parametrico
cubrimiento_cuantil_samp <- c(0,0,0,0,0,0)
cubrimiento_normal_samp <- c(0,0,0,0,0,0)

longitud_cuantil_samp <- c(0,0,0,0,0,0)
longitud_normal_samp <- c(0,0,0,0,0,0)

estimadores <- c(0,0,0,0,0,0)

titas.boot.samp <- matrix(0 , nrow = 6, ncol = 2500)

N_2 <- 1000

for (k in 1:N_2){
  for (j in 1:6){
    muestra <- rbinom(n[j],1,p_func(SE0,SP0,THETA0))
    estimadores[j] <- (mean(muestra)+SP0-1)/(SE0+SP0-1)
    muestras_boot <- replicate(2500, sample(muestra, replace = TRUE))
    titas.boot.samp[j, ] <- (colMeans(muestras_boot) + SP0 - 1) / (SE0 + SP0 - 1)
  }
  for (h in 1:6){
    se <- sqrt(mean((titas.boot.samp[h,] - mean(titas.boot.samp[h,]))**2))
    intervalo_normal_samp <- c(estimadores[h]-qnorm(0.975)*se, estimadores[h]+qnorm(0.975)*se)
    longitud_normal_samp[h] <- longitud_normal_samp[h]+intervalo_normal_samp[2]-intervalo_normal_samp[1]
    cubrimiento_normal_samp[h] <- cubrimiento_normal_samp[h]+(0.25>=intervalo_normal_samp[1] && 0.25<=intervalo_normal_samp[2])
    
    titas.boot.samp[h,] <- sort(titas.boot.samp[h,])
    intervalo_cuantil_samp <- c(quantile(titas.boot.samp[h,], 0.025), quantile(titas.boot.samp[h,], 0.975))
    longitud_cuantil_samp[h] <- longitud_cuantil_samp[h]+intervalo_cuantil_samp[2]-intervalo_cuantil_samp[1]
    cubrimiento_cuantil_samp[h] <- cubrimiento_cuantil_samp[h]+(0.25>=intervalo_cuantil_samp[1] && 0.25<=intervalo_cuantil_samp[2])
  }
}

longitud_normal_samp <- longitud_normal_samp/N_2
cubrimiento_normal_samp <- cubrimiento_normal_samp/N_2
longitud_cuantil_samp <- longitud_cuantil_samp/N_2
cubrimiento_cuantil_samp <- cubrimiento_cuantil_samp/N_2

resultados_comparacion_samp <- data.frame(
  n = n,
  Cubrimiento_Percentil = cubrimiento_cuantil_samp,
  Longitud_Percentil = longitud_cuantil_samp,
  Cubrimiento_Normal = cubrimiento_normal_samp,
  Longitud_Normal = longitud_normal_samp
)

print(resultados_comparacion_samp)

# Metodo parametrico

cubrimiento_cuantil_param <- c(0,0,0,0,0,0)
cubrimiento_normal_param <- c(0,0,0,0,0,0)

longitud_cuantil_param <- c(0,0,0,0,0,0)
longitud_normal_param <- c(0,0,0,0,0,0)

estimadores_param <- c(0,0,0,0,0,0)

titas.boot.param <- matrix(0 , nrow = 6, ncol = 2500)

for (k in 1:N_2){
  for (j in 1:6){
    muestra <- rbinom(n[j],1,p_func(SE0,SP0,THETA0))
    estimadores[j] <- (mean(muestra)+SP0-1)/(SE0+SP0-1)
    muestras_boot <- replicate(2500, sample(muestra, replace = TRUE))
    titas.boot.param[j, ] <- (colMeans(muestras_boot) + SP0 - 1) / (SE0 + SP0 - 1)
  }
  for (h in 1:6){
    se <- sqrt(mean((titas.boot.param[h,] - mean(titas.boot.param[h,]))**2))
    intervalo_normal_param <- c(estimadores_param[h]-qnorm(0.975)*se, estimadores_param[h]+qnorm(0.975)*se)
    longitud_normal_param[h] <- longitud_normal_param[h]+intervalo_normal_param[2]-intervalo_normal_param[1]
    cubrimiento_normal_parm[h] <- cubrimiento_normal_param[h]+(0.25>=intervalo_normal_param[1] && 0.25<=intervalo_normal_param[2])
    
    titas.boot.param[h,]<-sort(titas.boot.param[h,])
    intervalo_cuantil_param<-c(quantile(titas.boot.param[h,], 0.025), quantile(titas.boot.param[h,], 0.975))
    longitud_cuantil_param[h]<-longitud_cuantil_param[h]+intervalo_cuantil_param[2]-intervalo_cuantil_param[1]
    cubrimiento_cuantil_param[h]<-cubrimiento_cuantil_param[h]+(0.25>=intervalo_cuantil_param[1] && 0.25<=intervalo_cuantil_param[2])
  }
}

longitud_cuantil_param <- longitud_cuantil_param/N_2
cubrimiento_cuantil_param <- cubrimiento_cuantil_param/N_2
longitud_normal_param <- longitud_normal_param/N_2
cubrimiento_normal_param <- cubrimiento_normal_param/N_2

resultados_comparacion_param <- data.frame(
  n = n,
  Cubrimiento_Percentil = cubrimiento_cuantil_param,
  Longitud_Percentil = longitud_cuantil_param,
  Cubrimiento_Normal = cubrimiento_normal_param,
  Longitud_Normal = longitud_normal_param
)

print(resultados_comparacion_param)



# ----- Ejercicio 2.3: estimador truncado -----

n_trunc <- c(10,100,1000)

est_trunc <- function(estimador) {
  if (estimador>1){
    return(1)
  }
  if (estimador<0){
    return(0)
  }
  else {
    return(estimador)
  }
}

titas.boot.samp.trunc <- matrix(0 , nrow = 3, ncol = 2500)
estimadores.samp <- c(0,0,0)
titas.boot.param.trunc <- matrix(0 , nrow = 3, ncol = 2500)
estimadores.param <- c(0,0,0)

for (k in 1:N_2){
  for (j in 1:3){
    muestra <- rbinom(n_trunc[j],1,p_func(SE0,SP0,THETA0))
    estimadores.samp[j] <- (mean(muestra)+SP0-1)/(SE0+SP0-1)
    for (i in 1:2500){
      muestra_boot <- sample(muestra, replace=TRUE)
      theta_hat <- (mean(muestra_boot)+SP0-1)/(SE0+SP0-1)
      titas.boot.samp.trunc[j,i] <- theta_hat
    }
  }
}

for (k in 1:N_2){
  for (j in 1:3){
    muestra <- rbinom(n_trunc[j],1,p_func(SE0,SP0,THETA0))
    estimadores.param[j] <- (mean(muestra)+SP0-1)/(SE0+SP0-1)
    for (i in 1:2500){
      muestra_boot <- sample(muestra, replace=TRUE)
      theta_hat <- (mean(muestra_boot)+SP0-1)/(SE0+SP0-1)
      titas.boot.param.trunc[j,i] <- theta_hat
    }
  }
}
