rm(list = ls())
cat("\014")

# ============================================================
# LIBRERIAS
# ============================================================
library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(tidyr)
library(sandwich)   # errores robustos HC3
library(lmtest)     # coeftest con errores robustos

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_panel    <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")
ruta_resumen  <- file.path(RUTA_BASE, "output", "segunda_etapa_resumen.txt")
ruta_graficas <- file.path(RUTA_BASE, "output", "graficas")

dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# PARAMETROS
# ============================================================
REZAGO_VAR      <- 3     # mismo p que en pvar_estimacion.R
HORIZONTE_IRF   <- 12    # meses de acumulacion para el escalar de segunda etapa
SEMILLA         <- 2024

COLOR_PRINCIPAL <- "#2C5F8A"
COLOR_ACENTO    <- "#C94040"
COLOR_NEUTRO    <- "#A8A8A8"

# ============================================================
# LECTURA
# ============================================================
cat("Cargando panel completo...\n")

panel <- read_csv(ruta_panel, show_col_types = FALSE,
                  col_types = cols(cve_ent = col_integer())) %>%
  mutate(fecha = as_date(fecha)) %>%
  arrange(cve_ent, periodo)

cat("Panel:", nrow(panel), "obs |",
    n_distinct(panel$cve_ent), "estados |",
    n_distinct(panel$periodo), "periodos\n\n")

# ============================================================
# PARTIALLING OUT (replica exacta de pvar_estimacion.R)
# ============================================================
nombres_dummies <- paste0("d_mes", str_pad(2:12, 2, pad = "0"))
exogenas        <- c("d_tc", "cetes28", "expectativa", "d_wti", nombres_dummies)

proyectar_serie <- function(y, datos_estado, exogenas) {
  residuos <- rep(NA_real_, nrow(datos_estado))
  filas_ok <- complete.cases(datos_estado[, exogenas]) & !is.na(y)
  if (sum(filas_ok) > length(exogenas) + 1) {
    X    <- as.matrix(cbind(1, datos_estado[filas_ok, exogenas]))
    y_ok <- y[filas_ok]
    coef <- tryCatch(
      solve(t(X) %*% X) %*% t(X) %*% y_ok,
      error = function(e) NULL
    )
    if (!is.null(coef)) {
      residuos[filas_ok] <- y_ok - X %*% coef
    }
  }
  residuos
}

cat("Aplicando partialling out de exogenas...\n")

panel_limpio <- panel %>%
  arrange(cve_ent, periodo) %>%
  group_by(cve_ent) %>%
  group_modify(function(d, key) {
    d$h_r  <- proyectar_serie(d$h,  d, exogenas)
    d$g_r  <- proyectar_serie(d$g,  d, exogenas)
    d$pi_r <- proyectar_serie(d$pi, d, exogenas)
    d
  }) %>%
  ungroup()

# ============================================================
# FUNCIONES AUXILIARES PARA VAR INDIVIDUAL POR OLS
# ============================================================

# Construye la matriz de regresores con rezagos de las tres variables.
# Devuelve una lista con la matriz X y el vector de indices validos.
construir_rezagos <- function(datos_estado, p) {
  y_mat  <- as.matrix(datos_estado[, c("h_r", "g_r", "pi_r")])
  n_obs  <- nrow(y_mat)
  n_vars <- 3

  # Filas para las que se pueden construir todos los p rezagos
  primer_t <- p + 1
  n_validos <- n_obs - p

  if (n_validos < n_vars * p + 10) return(NULL)

  Y <- y_mat[(p + 1):n_obs, ]   # variable dependiente (n_validos x 3)

  # Construir X: intercepto + p bloques de n_vars columnas cada uno
  X_cols <- vector("list", p)
  for (lag in 1:p) {
    X_cols[[lag]] <- y_mat[(p + 1 - lag):(n_obs - lag), ]
  }
  X <- cbind(1, do.call(cbind, X_cols))

  # Eliminar filas con NA en Y o X
  filas_ok <- complete.cases(Y) & complete.cases(X)
  if (sum(filas_ok) < ncol(X) + 10) return(NULL)

  list(Y = Y[filas_ok, ], X = X[filas_ok, ], n_usadas = sum(filas_ok))
}

# Estima el VAR por OLS ecuacion por ecuacion y devuelve la
# matriz de coeficientes (filas = ecuaciones, columnas = regresores).
estimar_var_ols <- function(datos_estado, p) {
  partes <- construir_rezagos(datos_estado, p)
  if (is.null(partes)) return(NULL)

  Y <- partes$Y
  X <- partes$X

  # OLS para cada ecuacion por separado
  coef_mat <- matrix(NA_real_, nrow = 3, ncol = ncol(X))
  for (j in 1:3) {
    coef_j <- tryCatch(
      as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y[, j]),
      error = function(e) NULL
    )
    if (!is.null(coef_j)) coef_mat[j, ] <- coef_j
  }
  rownames(coef_mat) <- c("h_r", "g_r", "pi_r")
  coef_mat
}

# Construye el companion matrix de un VAR(p) con n_vars variables.
# La matriz de coeficientes tiene forma n_vars x (1 + n_vars*p).
# La companion matrix tiene dimension (n_vars*p) x (n_vars*p).
companion_matrix <- function(coef_mat, p) {
  n_vars <- nrow(coef_mat)
  # Coeficientes sin intercepto, organizados en p bloques de n_vars columnas
  A_cols <- coef_mat[, 2:ncol(coef_mat)]   # n_vars x (n_vars*p)

  # Extraer A_1, ..., A_p (cada una n_vars x n_vars)
  A_list <- vector("list", p)
  for (lag in 1:p) {
    idx_ini <- (lag - 1) * n_vars + 1
    idx_fin <- lag * n_vars
    A_list[[lag]] <- A_cols[, idx_ini:idx_fin]
  }

  # Construir companion matrix en bloques
  dim_comp <- n_vars * p
  C <- matrix(0, nrow = dim_comp, ncol = dim_comp)
  C[1:n_vars, ] <- A_cols

  if (p > 1) {
    C[(n_vars + 1):dim_comp, 1:(n_vars * (p - 1))] <-
      diag(n_vars * (p - 1))
  }
  C
}

# Calcula la IRF ortogonalizada de pi ante un choque en h usando
# la descomposicion de Cholesky con ordenamiento h -> g -> pi.
# Devuelve una matriz de (horizonte + 1) filas x n_vars columnas.
calcular_irf_ols <- function(coef_mat, sigma_mat, p, horizonte,
                              n_vars = 3) {
  # Descomposicion de Cholesky de la matriz de covarianza de residuos
  P_chol <- tryCatch(t(chol(sigma_mat)), error = function(e) NULL)
  if (is.null(P_chol)) return(NULL)

  C <- companion_matrix(coef_mat, p)
  dim_comp <- nrow(C)

  # IRF ortogonalizada: Phi_s = C^s restringido a las primeras n_vars filas
  # El impacto contemporaneo es P_chol
  irf_mat <- matrix(NA_real_, nrow = horizonte + 1, ncol = n_vars)

  C_potencia <- diag(dim_comp)   # C^0 = identidad
  for (s in 0:horizonte) {
    Phi_s <- C_potencia[1:n_vars, 1:n_vars]
    irf_mat[s + 1, ] <- Phi_s %*% P_chol[, 1]   # columna 1 = choque en h
    C_potencia <- C_potencia %*% C
  }
  irf_mat
}

# Estima la matriz de covarianza de los residuos de la regresion.
calcular_sigma <- function(datos_estado, coef_mat, p) {
  partes <- construir_rezagos(datos_estado, p)
  if (is.null(partes)) return(NULL)

  Y <- partes$Y
  X <- partes$X

  U <- Y - X %*% t(coef_mat)   # matriz de residuos (n_obs x n_vars)

  n_obs    <- nrow(U)
  n_params <- ncol(coef_mat)

  # Corriger los grados de libertad
  sigma <- (t(U) %*% U) / (n_obs - n_params)
  sigma
}

# ============================================================
# ETAPA 1: IRF INDIVIDUAL POR ESTADO
# ============================================================
cat("============================================================\n")
cat("ETAPA 1: VAR INDIVIDUAL POR ESTADO\n")
cat("============================================================\n\n")

estados_vec <- sort(unique(panel_limpio$cve_ent))
n_estados   <- length(estados_vec)

# Almacena la IRF de pi ante h para cada estado (horizonte + 1 puntos)
irf_por_estado    <- vector("list", n_estados)
names(irf_por_estado) <- as.character(estados_vec)

# Almacena el escalar: IRF acumulada de pi ante h a HORIZONTE_IRF meses
escalar_irf       <- rep(NA_real_, n_estados)
modulo_max        <- rep(NA_real_, n_estados)

for (i in seq_along(estados_vec)) {
  cve <- estados_vec[i]
  datos_estado <- panel_limpio %>%
    filter(cve_ent == cve) %>%
    arrange(periodo) %>%
    as.data.frame()

  coef_mat <- tryCatch(
    estimar_var_ols(datos_estado, REZAGO_VAR),
    error = function(e) NULL
  )

  if (is.null(coef_mat) || any(is.na(coef_mat))) {
    cat(sprintf("  Estado %02d: estimacion fallida\n", cve))
    next
  }

  sigma_mat <- tryCatch(
    calcular_sigma(datos_estado, coef_mat, REZAGO_VAR),
    error = function(e) NULL
  )

  if (is.null(sigma_mat) || any(is.na(sigma_mat))) {
    cat(sprintf("  Estado %02d: sigma fallida\n", cve))
    next
  }

  # Verificar estabilidad antes de calcular IRF
  C <- tryCatch(companion_matrix(coef_mat, REZAGO_VAR), error = function(e) NULL)
  if (!is.null(C)) {
    vp <- tryCatch(Mod(eigen(C, only.values = TRUE)$values),
                   error = function(e) rep(NA_real_, 1))
    modulo_max[i] <- max(vp, na.rm = TRUE)
    if (modulo_max[i] >= 1.05) {
      cat(sprintf("  Estado %02d: modulo maximo = %.3f (inestable, IRF no calculada)\n",
                  cve, modulo_max[i]))
      next
    }
  }

  irf_mat <- tryCatch(
    calcular_irf_ols(coef_mat, sigma_mat, REZAGO_VAR, HORIZONTE_IRF),
    error = function(e) NULL
  )

  if (is.null(irf_mat)) {
    cat(sprintf("  Estado %02d: IRF fallida\n", cve))
    next
  }

  # La columna 3 corresponde a pi (tercera variable en h -> g -> pi)
  irf_pi <- irf_mat[, 3]
  irf_por_estado[[as.character(cve)]] <- irf_pi

  # Respuesta acumulada de pi a HORIZONTE_IRF meses
  escalar_irf[i] <- sum(irf_pi[1:(HORIZONTE_IRF + 1)], na.rm = TRUE)

  cat(sprintf("  Estado %02d: IRF acumulada a %d meses = %.6f\n",
              cve, HORIZONTE_IRF, escalar_irf[i]))
}

cat("\nEstados con IRF calculada:",
    sum(!is.na(escalar_irf)), "de", n_estados, "\n\n")

# ============================================================
# GRAFICA g18: IRF DE PI ANTE H POR ESTADO (banda de cuantiles)
# ============================================================

# Construir data frame largo con las IRFs por estado
irf_lista_validas <- Filter(Negate(is.null), irf_por_estado)

if (length(irf_lista_validas) > 0) {
  irf_df_lista <- lapply(names(irf_lista_validas), function(cve) {
    data.frame(
      cve_ent   = as.integer(cve),
      horizonte = 0:HORIZONTE_IRF,
      irf_pi    = irf_lista_validas[[cve]]
    )
  })
  irf_todos <- do.call(rbind, irf_df_lista)

  cuantiles_irf <- irf_todos %>%
    group_by(horizonte) %>%
    summarise(
      mediana = median(irf_pi, na.rm = TRUE),
      p25     = quantile(irf_pi, 0.25, na.rm = TRUE),
      p75     = quantile(irf_pi, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  g18 <- ggplot(cuantiles_irf, aes(x = horizonte)) +
    geom_ribbon(aes(ymin = p25, ymax = p75),
                fill = COLOR_PRINCIPAL, alpha = 0.20) +
    geom_line(aes(y = mediana), color = COLOR_PRINCIPAL, linewidth = 0.9) +
    geom_hline(yintercept = 0, color = COLOR_NEUTRO,
               linewidth = 0.5, linetype = "dashed") +
    scale_x_continuous(breaks = 0:HORIZONTE_IRF) +
    labs(
      title    = "Respuesta de la inflacion ante un choque de violencia: distribucion entre estados",
      subtitle = paste0("VAR(", REZAGO_VAR, ") individual por estado. ",
                        "Banda: percentiles 25 y 75. Linea: mediana. ",
                        "Horizonte ", HORIZONTE_IRF, " meses"),
      x = "Meses despues del choque",
      y = "Respuesta de pi"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(ruta_graficas, "g18_segunda_etapa_irf_estados.png"),
         g18, width = 9, height = 5, dpi = 150)
  cat("Grafica g18 guardada.\n")
}

# ============================================================
# CARACTERISTICAS ESTRUCTURALES PARA LA SEGUNDA ETAPA
#
# Para obtener resultados comparables a Verdugo-Yepes et al. se
# usan dos caracteristicas como punto de partida:
#
#   1. informalidad_promedio: tasa promedio de ocupacion en el
#      sector informal por estado (ENOE 2005-2024). Si no se
#      dispone de estos datos, se usa el promedio nacional como
#      constante y la columna aparece como NA en la regresion.
#
#   2. violencia_promedio: tasa promedio de narco-homicidios
#      por 100 000 habitantes durante 2004-2024, construida
#      directamente del panel ya disponible.
#
# Para agregar regressores propios, construir un data frame
# adicional con columna cve_ent (integer) y las columnas
# deseadas, luego unirlo al data frame "caract" via left_join.
# ============================================================

# Caracteristica 2: violencia promedio (disponible en el panel)
violencia_promedio <- panel_limpio %>%
  group_by(cve_ent) %>%
  summarise(violencia_promedio = mean(h, na.rm = TRUE), .groups = "drop")

# Caracteristica 1: informalidad promedio por estado.
# Se carga desde un CSV externo si esta disponible en data/raw/.
# El archivo esperado se llama informalidad_estados.csv y debe
# tener columnas cve_ent (integer) e informalidad (proporcion 0-1).
# Si no existe, la columna se llena con NA y la regresion omitira
# ese regresor o el usuario puede completarla despues.
ruta_informalidad <- file.path(RUTA_BASE, "data", "raw",
                               "informalidad_estados.csv")

if (file.exists(ruta_informalidad)) {
  informalidad_df <- read_csv(ruta_informalidad,
                              col_types = cols(cve_ent = col_integer()),
                              show_col_types = FALSE)
  cat("Informalidad cargada desde:", ruta_informalidad, "\n")
} else {
  cat("Archivo informalidad_estados.csv no encontrado.\n")
  cat("La columna informalidad_promedio se inicializa como NA.\n")
  cat("Para completar: crear data/raw/informalidad_estados.csv\n")
  cat("  con columnas cve_ent (integer) e informalidad (0-1).\n\n")
  informalidad_df <- data.frame(
    cve_ent              = estados_vec,
    informalidad_promedio = NA_real_
  )
}

# Unir caracteristicas con los escalares de IRF
tabla_segunda_etapa <- data.frame(
  cve_ent   = estados_vec,
  irf_acum  = escalar_irf
) %>%
  left_join(violencia_promedio, by = "cve_ent") %>%
  left_join(
    informalidad_df %>%
      rename(informalidad_promedio = informalidad_promedio),
    by = "cve_ent"
  ) %>%
  filter(!is.na(irf_acum))

cat("\nTabla segunda etapa (primeras filas):\n")
print(head(tabla_segunda_etapa))
cat("Observaciones validas:", nrow(tabla_segunda_etapa), "\n\n")

# ============================================================
# ETAPA 2: REGRESION MCO CON ERRORES ROBUSTOS
# ============================================================
cat("============================================================\n")
cat("ETAPA 2: REGRESION TRANSVERSAL MCO\n")
cat("============================================================\n\n")

# Regresion sobre violencia promedio (siempre disponible)
reg_violencia <- lm(irf_acum ~ violencia_promedio,
                    data = tabla_segunda_etapa)
cat("Regresion: IRF acumulada ~ violencia_promedio\n")
print(coeftest(reg_violencia, vcov = vcovHC(reg_violencia, type = "HC3")))
cat("\n")

# Regresion sobre informalidad (si esta disponible)
tiene_informalidad <- !all(is.na(tabla_segunda_etapa$informalidad_promedio))
if (tiene_informalidad) {
  reg_informalidad <- lm(irf_acum ~ informalidad_promedio,
                         data = tabla_segunda_etapa)
  cat("Regresion: IRF acumulada ~ informalidad_promedio\n")
  print(coeftest(reg_informalidad, vcov = vcovHC(reg_informalidad, type = "HC3")))
  cat("\n")

  # Regresion multivariada
  reg_ambas <- lm(irf_acum ~ violencia_promedio + informalidad_promedio,
                  data = tabla_segunda_etapa)
  cat("Regresion: IRF acumulada ~ violencia_promedio + informalidad_promedio\n")
  print(coeftest(reg_ambas, vcov = vcovHC(reg_ambas, type = "HC3")))
  cat("\n")
}

# ============================================================
# GRAFICA g19: SCATTER IRF ACUMULADA vs. VIOLENCIA PROMEDIO
# ============================================================
g19 <- ggplot(tabla_segunda_etapa,
              aes(x = violencia_promedio, y = irf_acum)) +
  geom_point(color = COLOR_PRINCIPAL, size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = COLOR_ACENTO,
              fill = COLOR_ACENTO, alpha = 0.15, linewidth = 0.8) +
  labs(
    title    = "Segunda etapa: respuesta acumulada de la inflacion y violencia media",
    subtitle = paste0("IRF acumulada de pi ante choque en h a ", HORIZONTE_IRF,
                      " meses vs. tasa media de narco-homicidios 2004-2024"),
    x = "Tasa media de narco-homicidios (log, por 100 000 hab.)",
    y = paste0("IRF acumulada de pi a ", HORIZONTE_IRF, " meses")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(ruta_graficas, "g19_segunda_etapa_scatter_violencia.png"),
       g19, width = 8, height = 5, dpi = 150)
cat("Grafica g19 guardada.\n")

# Grafica g20 solo si hay datos de informalidad
if (tiene_informalidad) {
  g20 <- ggplot(tabla_segunda_etapa,
                aes(x = informalidad_promedio, y = irf_acum)) +
    geom_point(color = COLOR_PRINCIPAL, size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = COLOR_ACENTO,
                fill = COLOR_ACENTO, alpha = 0.15, linewidth = 0.8) +
    labs(
      title    = "Segunda etapa: respuesta acumulada de la inflacion e informalidad",
      subtitle = paste0("IRF acumulada de pi ante choque en h a ", HORIZONTE_IRF,
                        " meses vs. tasa media de informalidad laboral (ENOE)"),
      x = "Tasa media de informalidad laboral (proporcion)",
      y = paste0("IRF acumulada de pi a ", HORIZONTE_IRF, " meses")
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(ruta_graficas, "g20_segunda_etapa_scatter_informalidad.png"),
         g20, width = 8, height = 5, dpi = 150)
  cat("Grafica g20 guardada.\n")
}

# ============================================================
# EXPORTACION DEL RESUMEN
# ============================================================
sink(ruta_resumen)

cat("============================================================\n")
cat("SEGUNDA ETAPA: VAR INDIVIDUAL + REGRESION TRANSVERSAL\n")
cat("Panel: 32 estados, enero 2004 - diciembre 2024\n")
cat("VAR(p) individual por OLS, p =", REZAGO_VAR, "\n")
cat("Ordenamiento Cholesky: h -> g -> pi\n")
cat("IRF acumulada de pi ante choque en h, horizonte:", HORIZONTE_IRF, "meses\n")
cat("Segunda etapa: MCO con errores robustos HC3\n")
cat("============================================================\n\n")

cat("Estados con IRF valida:", sum(!is.na(escalar_irf)), "\n\n")

cat("Tabla escalares IRF por estado:\n")
print(tabla_segunda_etapa)

cat("\n\nRegresion: IRF acumulada ~ violencia_promedio\n")
print(summary(reg_violencia))
print(coeftest(reg_violencia, vcov = vcovHC(reg_violencia, type = "HC3")))

if (tiene_informalidad) {
  cat("\n\nRegresion: IRF acumulada ~ informalidad_promedio\n")
  print(summary(reg_informalidad))
  print(coeftest(reg_informalidad, vcov = vcovHC(reg_informalidad, type = "HC3")))

  cat("\n\nRegresion multivariada\n")
  print(summary(reg_ambas))
  print(coeftest(reg_ambas, vcov = vcovHC(reg_ambas, type = "HC3")))
}

sink()

cat("\nResumen exportado a:", ruta_resumen, "\n")
