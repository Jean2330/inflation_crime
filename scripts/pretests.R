rm(list = ls())
cat("\014")

# ============================================================
# LIBRERIAS
# ============================================================
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(plm)    # purtest(): prueba IPS de raiz unitaria en panel
library(urca)   # ca.po():   prueba de cointegracion de Pedroni

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_panel   <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")
ruta_resumen <- file.path(RUTA_BASE, "output", "pretests_resumen.txt")

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
# PREPARACION COMO PDATA.FRAME
# purtest() requiere un pdata.frame indexado por estado y periodo.
# Los identificadores deben ser character o factor, no integer.
# ============================================================
panel_plm <- pdata.frame(
  panel %>% mutate(cve_ent = as.character(cve_ent)),
  index = c("cve_ent", "periodo")
)

# ============================================================
# PRUEBA IPS DE RAIZ UNITARIA
# H0: todas las series son I(1)
# H1: al menos una serie es estacionaria
# Se prueba con intercepto y con intercepto mas tendencia.
# lags = "AIC" selecciona el orden del ADF individual por AIC.
# ============================================================
cat("============================================================\n")
cat("PRUEBA IPS DE RAIZ UNITARIA (plm::purtest)\n")
cat("H0: todas las series tienen raiz unitaria\n")
cat("============================================================\n\n")

# Funcion auxiliar que corre IPS con las dos especificaciones
# y devuelve una lista con ambos resultados.
correr_ips <- function(serie, nombre) {
  cat("--- Variable:", nombre, "---\n")
  
  ips_c <- purtest(serie, test = "ips", exo = "intercept",
                   lags = "AIC", pmax = 12)
  cat("Con intercepto:\n")
  print(summary(ips_c))
  
  ips_ct <- purtest(serie, test = "ips", exo = "trend",
                    lags = "AIC", pmax = 12)
  cat("Con intercepto y tendencia:\n")
  print(summary(ips_ct))
  
  cat("\n")
  list(intercepto = ips_c, tendencia = ips_ct)
}

resultados_ips <- list(
  h  = correr_ips(panel_plm$h,  "h  (violencia)"),
  g  = correr_ips(panel_plm$g,  "g  (brecha ITAEE)"),
  pi = correr_ips(panel_plm$pi, "pi (inflacion)")
)

# ============================================================
# IPS SOBRE PRIMERAS DIFERENCIAS
# Confirma el orden de integracion: si la serie en niveles no
# rechaza H0 pero la primera diferencia si lo hace, es I(1).
#
# Se usa plm::diff() directamente sobre el pdata.frame de niveles.
# Esto evita el problema de NAs y huecos en el indice: diff()
# sobre un pdata.frame respeta la estructura de panel y devuelve
# una pseries sin huecos, que purtest acepta sin error.
# ============================================================
cat("============================================================\n")
cat("IPS EN PRIMERAS DIFERENCIAS\n")
cat("============================================================\n\n")

panel_plm$dh  <- diff(panel_plm$h)
panel_plm$dg  <- diff(panel_plm$g)
panel_plm$dpi <- diff(panel_plm$pi)

resultados_ips_dif <- list(
  dh  = correr_ips(panel_plm$dh,  "dh  (primera diferencia de h)"),
  dg  = correr_ips(panel_plm$dg,  "dg  (primera diferencia de g)"),
  dpi = correr_ips(panel_plm$dpi, "dpi (primera diferencia de pi)")
)

# ============================================================
# PRUEBA DE COINTEGRACION DE PEDRONI
# ca.po() de urca implementa Phillips y Ouliaris (1990), que es
# la base sobre la que Pedroni (1999, 2004) extiende al panel.
# Se estima la regresion de cointegracion estado por estado y
# se acumulan los estadisticos individuales para el test de panel.
#
# Estrategia: para cada estado, correr ca.po() con pi como
# variable dependiente y (h, g) como regresoras. Se reportan
# los estadisticos individuales y el promedio del panel.
#
# H0: no existe cointegracion
# H1: existe cointegracion
# ============================================================
cat("============================================================\n")
cat("COINTEGRACION DE PEDRONI ESTADO POR ESTADO (urca::ca.po)\n")
cat("H0: no existe cointegracion entre h, g y pi\n")
cat("============================================================\n\n")

lista_estados <- sort(unique(panel$cve_ent))

estadisticos_po <- data.frame(
  cve_ent = integer(0),
  estadistico = numeric(0),
  valor_critico_5pct = numeric(0)
)

for (estado in lista_estados) {
  datos_estado <- panel %>%
    filter(cve_ent == estado) %>%
    arrange(periodo) %>%
    select(pi, h, g) %>%
    drop_na()
  
  if (nrow(datos_estado) < 20) next
  
  resultado_po <- tryCatch(
    ca.po(as.matrix(datos_estado), type = "Pu", demean = "const"),
    error = function(e) NULL
  )
  
  if (!is.null(resultado_po)) {
    estad <- resultado_po@teststat
    crit  <- resultado_po@cval[2]   # valor critico al 5%
    estadisticos_po <- rbind(
      estadisticos_po,
      data.frame(cve_ent           = estado,
                 estadistico        = as.numeric(estad),
                 valor_critico_5pct = as.numeric(crit))
    )
  }
}

cat("Estadisticos de Phillips-Ouliaris por estado:\n")
print(estadisticos_po)

n_rechazan <- sum(estadisticos_po$estadistico >
                    estadisticos_po$valor_critico_5pct, na.rm = TRUE)
cat("\nEstados que rechazan H0 al 5%:", n_rechazan,
    "de", nrow(estadisticos_po), "\n\n")

# Promedio de panel segun la logica de Pedroni: si la mayoria
# de los estados rechazan H0 individualmente, hay evidencia de
# cointegracion en el panel.
cat("Proporcion de estados con cointegracion:",
    round(n_rechazan / nrow(estadisticos_po), 3), "\n\n")

# ============================================================
# EXPORTACION DEL RESUMEN
# ============================================================
sink(ruta_resumen)

cat("============================================================\n")
cat("RESUMEN DE PRUEBAS DE RAIZ UNITARIA Y COINTEGRACION\n")
cat("Panel: 32 estados, enero 2004 - diciembre 2024\n")
cat("IPS: Im, Pesaran y Shin (2003) via plm::purtest\n")
cat("Cointegracion: Phillips-Ouliaris via urca::ca.po (logica Pedroni)\n")
cat("============================================================\n\n")

cat("Criterio de interpretacion:\n")
cat("IPS rechaza H0 -> la variable es estacionaria en el panel.\n")
cat("Si alguna es I(1) y la mayoria de estados muestran cointegracion,\n")
cat("la estimacion en niveles sigue siendo valida.\n\n")

for (nombre in names(resultados_ips)) {
  cat("Variable:", nombre, "\n")
  cat("IPS con intercepto:\n")
  print(summary(resultados_ips[[nombre]]$intercepto))
  cat("IPS con tendencia:\n")
  print(summary(resultados_ips[[nombre]]$tendencia))
  cat("\n")
}

cat("IPS en primeras diferencias:\n")
for (nombre in names(resultados_ips_dif)) {
  cat("Variable:", nombre, "\n")
  print(summary(resultados_ips_dif[[nombre]]$intercepto))
  cat("\n")
}

cat("Phillips-Ouliaris por estado (pi ~ h + g):\n")
print(estadisticos_po)
cat("\nEstados que rechazan H0 al 5%:", n_rechazan,
    "de", nrow(estadisticos_po), "\n")
cat("Proporcion:", round(n_rechazan / nrow(estadisticos_po), 3), "\n")

sink()

cat("Resumen exportado a:", ruta_resumen, "\n")