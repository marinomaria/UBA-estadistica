library(latex2exp) # Para gráficos

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

se_real = 0.9
sp_real = 0.95
theta_real = 0.25

# Definimos la función p
p <- function(se, sp, theta) {
  return(theta * (se + sp - 1) + 1 - sp)
}

par(mfrow = c(1, 3)) 

curve(p(x, sp_real, theta_real), 
      from = 0, to = 1, 
      main = TeX("Variación en $Se$"),
      xlab = TeX("$Se$"), ylab = TeX("$p(Se, Sp, \theta)$"),
      col = "blue", lwd = 2)
grid()

curve(p(se_real, x, theta_real), 
      from = 0, to = 1, 
      main = TeX("Variación en $Sp$"),
      xlab = TeX("$Sp$"), ylab = TeX("$p(Se, Sp, \theta)$"),
      col = "green", lwd = 2)
grid()

curve(p(se_real, sp_real, x), 
      from = 0, to = 1, 
      main = TeX("Variación en $\theta$"),
      xlab = TeX("$\theta$"), ylab = TeX("$p(Se, Sp, \theta)$"),
      col = "purple", lwd = 2)
grid()




