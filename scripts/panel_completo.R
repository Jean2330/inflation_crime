rm(list = ls())
cat("\014")

# ============================================================
# LIBRERIAS
# ============================================================
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_inpc  <- file.path(RUTA_BASE, "data", "processed", "inpc_ciudades.csv")
ruta_itaee <- file.path(RUTA_BASE, "data", "processed", "itaee.csv")
ruta_narco <- file.path(RUTA_BASE, "data", "raw",       "narco_homicidios.csv")
ruta_macro <- file.path(RUTA_BASE, "data", "processed", "macro_controles.csv")
ruta_panel <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")

# ============================================================
# PARAMETROS DEL FILTRO DE HAMILTON
# Hamilton (2018): regresion de y(t) sobre y(t-h), ..., y(t-h-p+1)
# h = horizonte de proyeccion, p = numero de rezagos
# ============================================================
HAMILTON_H <- 8    # horizonte en trimestres equivalentes (8 meses aqui, frecuencia mensual)
HAMILTON_P <- 4    # numero de rezagos en la regresion auxiliar

# ============================================================
# FUNCION: FILTRO DE HAMILTON
# Calcula el ciclo como el residuo de una regresion OLS de
# y(t+h) sobre y(t), y(t-1), ..., y(t-p+1).
# Devuelve el ciclo alineado con el vector original.
# ============================================================
aplicar_hamilton <- function(y, h = HAMILTON_H, p = HAMILTON_P) {
  n <- length(y)
  ciclo <- rep(NA_real_, n)

  # El primer indice para el que la regresion es posible es h + p
  primer_t <- h + p

  if (n < primer_t + 1) return(ciclo)

  for (t in primer_t:n) {
    # Variable dependiente: y adelantado h periodos
    y_adelantado <- y[t]
    if (is.na(y_adelantado)) next

    # Variables explicativas: y(t-h), y(t-h-1), ..., y(t-h-p+1)
    inicio_rezago <- t - h
    fin_rezago    <- t - h - p + 1

    if (fin_rezago < 1) next

    rezagos <- y[inicio_rezago:fin_rezago]
    if (any(is.na(rezagos))) next

    # Regresion OLS con intercepto
    X <- cbind(1, matrix(rezagos, nrow = 1))
    # Se acumula el sistema de ecuaciones y se resuelve al final
    # (aqui se calcula residuo punto a punto para no requerir matrices grandes)
    # Esta implementacion es correcta para series de hasta ~300 obs.
    ciclo[t] <- NA_real_  # placeholder; se sobreescribe abajo
  }

  # Implementacion vectorizada correcta del filtro de Hamilton:
  # Para cada estado, se estima una sola regresion con todos los t validos.
  indices_y  <- (primer_t):n
  y_dep      <- y[indices_y]

  # Construir la matriz de regresores para todos los t a la vez
  mat_X <- matrix(NA_real_, nrow = length(indices_y), ncol = p + 1)
  mat_X[, 1] <- 1  # intercepto

  for (k in 1:p) {
    rezago_k <- indices_y - h - (k - 1)
    valores_k <- ifelse(rezago_k >= 1, y[pmax(rezago_k, 1)], NA_real_)
    valores_k[rezago_k < 1] <- NA_real_
    mat_X[, k + 1] <- valores_k
  }

  filas_completas <- complete.cases(mat_X, y_dep)

  if (sum(filas_completas) < p + 2) return(ciclo)

  X_ok <- mat_X[filas_completas, , drop = FALSE]
  y_ok <- y_dep[filas_completas]

  coef_hamilton <- solve(t(X_ok) %*% X_ok) %*% t(X_ok) %*% y_ok
  residuos      <- y_ok - X_ok %*% coef_hamilton

  ciclo[indices_y[filas_completas]] <- as.numeric(residuos)

  return(ciclo)
}

# ============================================================
# LECTURA DE LOS CUATRO DATASETS
# ============================================================
cat("Cargando datos...\n")

inpc <- read_csv(ruta_inpc, show_col_types = FALSE,
                 col_types = cols(cve_ent = col_integer())) %>%
  mutate(fecha = as_date(fecha))

itaee <- read_csv(ruta_itaee, show_col_types = FALSE,
                  col_types = cols(cve_ent = col_integer())) %>%
  mutate(fecha = as_date(fecha))

narco_raw <- read_csv(ruta_narco, show_col_types = FALSE) %>%
  mutate(
    cve_ent = as.integer(str_sub(str_pad(cve_inegi, 5, pad = "0"), 1, 2)),
    fecha   = dmy(fecha)
  )

macro <- read_csv(ruta_macro, show_col_types = FALSE) %>%
  mutate(fecha = as_date(fecha))

cat("Lectura completa.\n\n")

# ============================================================
# CONSTRUCCION DE LA VARIABLE DE VIOLENCIA (h)
# Se agrega narco_homicidios.csv al nivel estado-mes,
# se calcula la tasa por 100,000 habitantes y se aplica
# la transformacion logaritmica log(1 + tasa).
# ============================================================
narco_estatal <- narco_raw %>%
  group_by(cve_ent, fecha) %>%
  summarise(
    homicidios = sum(homicidios_narco, na.rm = TRUE),
    poblacion  = sum(poblacion_total,  na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    tasa_narco = homicidios / poblacion * 100000,
    h          = log(1 + tasa_narco)
  ) %>%
  select(cve_ent, fecha, h, tasa_narco)

# ============================================================
# CONSTRUCCION DE LA INFLACION MENSUAL (pi)
# pi = delta log(INPC): primera diferencia del logaritmo del nivel.
# Se ordena por estado y fecha antes de diferenciar.
# ============================================================
inpc_pi <- inpc %>%
  arrange(cve_ent, fecha) %>%
  group_by(cve_ent) %>%
  mutate(pi = log(inpc) - lag(log(inpc))) %>%
  ungroup() %>%
  select(cve_ent, entidad, fecha, inpc, pi)

# ============================================================
# FILTRO DE HAMILTON SOBRE EL LOG DEL ITAEE
# Se aplica estado por estado sobre el logaritmo del ITAEE.
# El ciclo resultante es la brecha del producto (g).
# ============================================================
cat("Aplicando filtro de Hamilton al ITAEE...\n")

itaee_con_brecha <- itaee %>%
  arrange(cve_ent, fecha) %>%
  group_by(cve_ent) %>%
  mutate(
    log_itaee = log(itaee),
    g         = aplicar_hamilton(log_itaee)
  ) %>%
  ungroup() %>%
  select(cve_ent, fecha, itaee, log_itaee, g)

cat("Filtro de Hamilton aplicado.\n\n")

# ============================================================
# CONSTRUCCION DE LAS VARIABLES EXOGENAS
# Se utilizan los controles macroeconomicos nacionales:
#   d_tc        = depreciacion cambiaria (log-diferencia TC FIX)
#   cetes28     = tasa de CETES a 28 dias (nivel)
#   expectativa = suma de las 12 expectativas mensuales de inflacion
#   d_wti       = cambio en precio del petroleo (log-diferencia WTI)
# ============================================================
macro_exogenas <- macro %>%
  arrange(fecha) %>%
  mutate(
    d_tc  = log(tc_fix)   - lag(log(tc_fix)),
    d_wti = log(petroleo) - lag(log(petroleo))
  ) %>%
  select(fecha, d_tc, cetes28, expectativa, d_wti)

# ============================================================
# DUMMIES DE MES PARA CONTROLAR ESTACIONALIDAD DEL INPC
# Se generan 11 dummies (mes 1 es la categoria base) para
# capturar patrones estacionales del nivel de precios.
# ============================================================
meses_secuencia <- sort(unique(inpc_pi$fecha))
dummies_mes <- data.frame(
  fecha = meses_secuencia,
  mes   = month(meses_secuencia)
)

for (m in 2:12) {
  nombre_dummy <- paste0("d_mes", str_pad(m, 2, pad = "0"))
  dummies_mes[[nombre_dummy]] <- as.integer(dummies_mes$mes == m)
}

dummies_mes <- select(dummies_mes, -mes)

# ============================================================
# UNION DE TODAS LAS FUENTES
# El panel resultante tiene una observacion por estado y mes.
# Se filtra al periodo 2004-2024.
# ============================================================
cat("Construyendo panel maestro...\n")

panel <- inpc_pi %>%
  inner_join(itaee_con_brecha, by = c("cve_ent", "fecha")) %>%
  left_join(narco_estatal,     by = c("cve_ent", "fecha")) %>%
  left_join(macro_exogenas,    by = "fecha") %>%
  left_join(dummies_mes,       by = "fecha") %>%
  filter(fecha >= as_date(FECHA_INICIO),
         fecha <= as_date(FECHA_FIN)) %>%
  arrange(cve_ent, fecha)

# ============================================================
# IDENTIFICADOR DE PERIODO (entero) REQUERIDO POR pvargmm
# pvargmm no acepta fechas: el tiempo debe ser un entero
# que indexa la posicion de cada mes en la secuencia del panel.
# ============================================================
fechas_ordenadas <- sort(unique(panel$fecha))
panel <- panel %>%
  mutate(periodo = match(fecha, fechas_ordenadas))

# ============================================================
# REORDENAMIENTO DE COLUMNAS
# cve_ent y periodo deben ser columnas 1 y 2 para pvargmm.
# ============================================================
panel <- panel %>%
  select(cve_ent, periodo, entidad, fecha,
         h, tasa_narco,
         g, log_itaee, itaee,
         pi, inpc,
         d_tc, cetes28, expectativa, d_wti,
         starts_with("d_mes"))

# ============================================================
# VERIFICACION DE BALANCE Y CALIDAD
# ============================================================
cat("Estados en el panel:", n_distinct(panel$cve_ent), "\n")
cat("Periodos en el panel:", n_distinct(panel$periodo), "\n")
cat("Observaciones totales:", nrow(panel), "\n")
cat("Esperado (32 x 252):", 32 * 252, "\n")

cat("\nNAs por variable clave:\n")
vars_clave <- c("h", "g", "pi", "d_tc", "cetes28", "expectativa", "d_wti")
for (v in vars_clave) {
  cat(" ", v, ":", sum(is.na(panel[[v]])), "\n")
}

cat("\nFechas:", format(min(panel$fecha)), "a", format(max(panel$fecha)), "\n\n")

# ============================================================
# EXPORTACION
# ============================================================
write_csv(panel, ruta_panel)
cat("Panel exportado a:", ruta_panel, "\n")
cat("Observaciones exportadas:", nrow(panel), "\n")
