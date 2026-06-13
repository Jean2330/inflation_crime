rm(list = ls())
cat("\014")

# ============================================================
# LIBRERIAS
# ============================================================
library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(panelvar)
library(ggplot2)
library(vars)
library(siebanxicor)

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_panel    <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")
ruta_resumen  <- file.path(RUTA_BASE, "output", "robustez_resumen.txt")
ruta_graficas <- file.path(RUTA_BASE, "output", "graficas")

dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# PARAMETROS
# Actualizar REZAGO_OPTIMO segun el resultado del MMSC en pvar_estimacion.R.
# ============================================================
REZAGO_OPTIMO    <- 3
HORIZONTE_IRF    <- 24
HORIZONTE_ETAPA2 <- 12     # meses para IRF acumulado en segunda etapa
SEMILLA          <- 2024
UMBRAL_VIOLENCIA <- 2.0    # tasa narco promedio por 100k a nivel estatal-mes

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
  arrange(cve_ent, periodo) %>%
  dplyr::select(cve_ent, periodo, everything()) %>%
  as.data.frame()

cat("Panel:", nrow(panel), "obs\n\n")

# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

nombres_dummies  <- paste0("d_mes", str_pad(2:12, 2, pad = "0"))
exogenas_base    <- c("d_tc", "cetes28", "expectativa", "d_wti", nombres_dummies)

# Proyecta una serie sobre las exogenas por OLS estado por estado.
# Intersecta filas sin NAs en exogenas y sin NAs en y para manejar
# los NAs iniciales de g producidos por el filtro de Hamilton.
proyectar_serie <- function(y, datos_estado, exogenas) {
  residuos <- rep(NA_real_, nrow(datos_estado))
  filas_ok <- complete.cases(datos_estado[, exogenas]) & !is.na(y)
  if (sum(filas_ok) > length(exogenas) + 1) {
    X    <- as.matrix(cbind(1, datos_estado[filas_ok, exogenas]))
    y_ok <- y[filas_ok]
    coef <- tryCatch(solve(t(X) %*% X) %*% t(X) %*% y_ok,
                     error = function(e) NULL)
    if (!is.null(coef)) residuos[filas_ok] <- y_ok - X %*% coef
  }
  residuos
}

# Aplica partialling out sobre un subpanel usando el vector de exogenas dado.
# Devuelve el panel con columnas h_r, g_r, pi_r agregadas.
aplicar_partialling <- function(datos, exogenas) {
  datos %>%
    arrange(cve_ent, periodo) %>%
    group_by(cve_ent) %>%
    group_modify(function(d, key) {
      d$h_r  <- proyectar_serie(d$h,  d, exogenas)
      d$g_r  <- proyectar_serie(d$g,  d, exogenas)
      d$pi_r <- proyectar_serie(d$pi, d, exogenas)
      d
    }) %>%
    ungroup() %>%
    as.data.frame()
}

# Aplica partialling out y estima PVAR-GMM con el vector endogeno dado.
estimar_pvar <- function(datos, endogenas, rezagos, exogenas = exogenas_base) {
  datos_limpios <- aplicar_partialling(datos, exogenas)
  tryCatch(
    pvargmm(
      dependent_vars   = endogenas,
      lags             = rezagos,
      transformation   = "fod",
      data             = datos_limpios,
      panel_identifier = c(1, 2),
      steps            = "onestep",
      collapse         = TRUE
    ),
    error = function(e) {
      cat("Error en estimacion:", conditionMessage(e), "\n")
      NULL
    }
  )
}

# Calcula IRFs ortogonalizadas con oirf().
calcular_irf <- function(modelo, horizonte, semilla) {
  set.seed(semilla)
  oirf(modelo, n.ahead = horizonte)
}

# Extrae y grafica la IRF de un par impulso-respuesta especifico.
# oirf() devuelve un objeto de clase pvaroirf donde las claves de nivel
# superior son las variables RESPUESTA. Cada elemento es una matriz de
# horizonte x n_vars columnas, donde las columnas son los choques en el
# orden del vector endogeno del modelo estimado.
graficar_irf_robustez <- function(irf_obj, endogenas_modelo,
                                  impulso, respuesta,
                                  titulo, subtitulo, ruta_salida) {
  col_impulso <- which(endogenas_modelo == impulso)
  
  if (length(col_impulso) == 0) {
    cat("Impulso no reconocido:", impulso, "\n")
    return(invisible(NULL))
  }
  if (!respuesta %in% names(irf_obj)) {
    cat("Respuesta no encontrada:", respuesta, "\n")
    cat("Disponibles:", paste(names(irf_obj), collapse = ", "), "\n")
    return(invisible(NULL))
  }
  
  mat <- irf_obj[[respuesta]]
  datos_irf <- data.frame(
    horizonte = 0:(nrow(mat) - 1),
    media     = mat[, col_impulso],
    lower     = mat[, col_impulso],
    upper     = mat[, col_impulso]
  )
  
  g <- ggplot(datos_irf, aes(x = horizonte)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = COLOR_PRINCIPAL, alpha = 0.2) +
    geom_line(aes(y = media), color = COLOR_PRINCIPAL, linewidth = 0.9) +
    geom_hline(yintercept = 0, color = COLOR_NEUTRO, linewidth = 0.5,
               linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, HORIZONTE_IRF, by = 6)) +
    labs(title = titulo, subtitle = subtitulo,
         x = "Meses despues del choque", y = "Respuesta") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(ruta_salida, g, width = 8, height = 5, dpi = 150)
  cat("Grafica guardada:", ruta_salida, "\n")
}

# Re-indexa periodo como entero consecutivo dentro del subpanel.
# pvargmm requiere que el indice de tiempo no tenga huecos.
reindexar_periodo <- function(datos) {
  fechas_unicas  <- sort(unique(datos$fecha))
  datos$periodo  <- match(datos$fecha, fechas_unicas)
  datos
}

# Prueba de causalidad de Granger manual para un par causa -> causada.
# Usa los errores estandar del primer paso (standard_error_first_step).
# W = sum((beta / se)^2) ~ Chi2(n_lags) bajo H0 de no causalidad.
probar_granger <- function(modelo, causa, causada) {
  vars_end  <- modelo$dependent_vars
  n_lags    <- modelo$lags
  n_vars    <- length(vars_end)
  col_causa  <- which(vars_end == causa)
  col_causada <- which(vars_end == causada)
  
  if (length(col_causa) == 0 || length(col_causada) == 0) {
    cat("Variable no encontrada:", causa, "->", causada, "\n")
    return(data.frame(causa = causa, causada = causada,
                      wald = NA_real_, df = n_lags, p_valor = NA_real_))
  }
  
  idx_coef <- seq(col_causa, n_vars * n_lags, by = n_vars)
  beta_sub <- modelo$coefficients[col_causada, idx_coef]
  se_sub   <- modelo$standard_error_first_step[col_causada, idx_coef]
  
  W_stat <- tryCatch(sum((beta_sub / se_sub)^2), error = function(e) NA_real_)
  p_val  <- 1 - pchisq(W_stat, df = n_lags)
  
  cat(sprintf("  %s -> %s: Wald = %.3f, df = %d, p = %.4f\n",
              causa, causada, W_stat, n_lags, p_val))
  data.frame(causa = causa, causada = causada,
             wald = W_stat, df = n_lags, p_valor = p_val)
}

endogenas_principal <- c("h_r", "g_r", "pi_r")
endogenas_irf       <- c("pi_r", "g_r", "h_r")   # reordenamiento solo para graficas

# ============================================================
# EJERCICIO 1: ESTADOS VIOLENTOS VS PACIFICOS
# ============================================================
cat("============================================================\n")
cat("EJERCICIO 1: ESTADOS VIOLENTOS VS PACIFICOS\n")
cat("============================================================\n\n")

tasa_promedio_estado <- panel %>%
  group_by(cve_ent) %>%
  summarise(tasa_media = mean(tasa_narco, na.rm = TRUE), .groups = "drop")

estados_violentos <- tasa_promedio_estado %>%
  filter(tasa_media >= UMBRAL_VIOLENCIA) %>% pull(cve_ent)

estados_pacificos <- tasa_promedio_estado %>%
  filter(tasa_media <  UMBRAL_VIOLENCIA) %>% pull(cve_ent)

cat("Estados violentos (tasa >=", UMBRAL_VIOLENCIA, "):", length(estados_violentos), "\n")
cat("Estados pacificos (tasa < ",  UMBRAL_VIOLENCIA, "):", length(estados_pacificos), "\n\n")

panel_violentos <- panel %>%
  filter(cve_ent %in% estados_violentos) %>%
  reindexar_periodo() %>% as.data.frame()

modelo_violentos <- estimar_pvar(panel_violentos, endogenas_principal, REZAGO_OPTIMO)

if (!is.null(modelo_violentos)) {
  modelo_violentos_irf <- estimar_pvar(panel_violentos, endogenas_irf, REZAGO_OPTIMO)
  irf_violentos <- calcular_irf(
    if (!is.null(modelo_violentos_irf)) modelo_violentos_irf else modelo_violentos,
    HORIZONTE_IRF, SEMILLA)
  graficar_irf_robustez(
    irf_violentos, endogenas_irf, "h_r", "pi_r",
    titulo      = "IRF en estados violentos: inflacion ante choque de violencia",
    subtitulo   = paste0("Estados con tasa narco >= ", UMBRAL_VIOLENCIA,
                         " por 100k. IRF ortogonalizada, p = ", REZAGO_OPTIMO),
    ruta_salida = file.path(ruta_graficas, "g13_irf_violentos.png")
  )
}

panel_pacificos <- panel %>%
  filter(cve_ent %in% estados_pacificos) %>%
  reindexar_periodo() %>% as.data.frame()

modelo_pacificos <- estimar_pvar(panel_pacificos, endogenas_principal, REZAGO_OPTIMO)

if (!is.null(modelo_pacificos)) {
  modelo_pacificos_irf <- estimar_pvar(panel_pacificos, endogenas_irf, REZAGO_OPTIMO)
  irf_pacificos <- calcular_irf(
    if (!is.null(modelo_pacificos_irf)) modelo_pacificos_irf else modelo_pacificos,
    HORIZONTE_IRF, SEMILLA)
  graficar_irf_robustez(
    irf_pacificos, endogenas_irf, "h_r", "pi_r",
    titulo      = "IRF en estados pacificos: inflacion ante choque de violencia",
    subtitulo   = paste0("Estados con tasa narco < ", UMBRAL_VIOLENCIA,
                         " por 100k. IRF ortogonalizada, p = ", REZAGO_OPTIMO),
    ruta_salida = file.path(ruta_graficas, "g14_irf_pacificos.png")
  )
}

# ============================================================
# EJERCICIO 2: SUB-PERIODOS
# ============================================================
cat("\n============================================================\n")
cat("EJERCICIO 2: SUB-PERIODOS\n")
cat("============================================================\n\n")

panel_sp1 <- panel %>%
  filter(fecha >= "2004-01-01", fecha <= "2012-12-01") %>%
  reindexar_periodo() %>% as.data.frame()

panel_sp2 <- panel %>%
  filter(fecha >= "2013-01-01", fecha <= "2024-12-01") %>%
  reindexar_periodo() %>% as.data.frame()

cat("Sub-periodo 1 (2004-2012):", nrow(panel_sp1), "obs\n")
cat("Sub-periodo 2 (2013-2024):", nrow(panel_sp2), "obs\n\n")

modelo_sp1 <- estimar_pvar(panel_sp1, endogenas_principal, REZAGO_OPTIMO)
if (!is.null(modelo_sp1)) {
  modelo_sp1_irf <- estimar_pvar(panel_sp1, endogenas_irf, REZAGO_OPTIMO)
  irf_sp1 <- calcular_irf(
    if (!is.null(modelo_sp1_irf)) modelo_sp1_irf else modelo_sp1,
    HORIZONTE_IRF, SEMILLA)
  graficar_irf_robustez(
    irf_sp1, endogenas_irf, "h_r", "pi_r",
    titulo      = "IRF 2004-2012: inflacion ante choque de violencia",
    subtitulo   = paste0("Primer sub-periodo. IRF ortogonalizada, p = ", REZAGO_OPTIMO),
    ruta_salida = file.path(ruta_graficas, "g15_irf_subperiodo1.png")
  )
}

modelo_sp2 <- estimar_pvar(panel_sp2, endogenas_principal, REZAGO_OPTIMO)
if (!is.null(modelo_sp2)) {
  modelo_sp2_irf <- estimar_pvar(panel_sp2, endogenas_irf, REZAGO_OPTIMO)
  irf_sp2 <- calcular_irf(
    if (!is.null(modelo_sp2_irf)) modelo_sp2_irf else modelo_sp2,
    HORIZONTE_IRF, SEMILLA)
  graficar_irf_robustez(
    irf_sp2, endogenas_irf, "h_r", "pi_r",
    titulo      = "IRF 2013-2024: inflacion ante choque de violencia",
    subtitulo   = paste0("Segundo sub-periodo. IRF ortogonalizada, p = ", REZAGO_OPTIMO),
    ruta_salida = file.path(ruta_graficas, "g16_irf_subperiodo2.png")
  )
}

# ============================================================
# EJERCICIO 3: RE-ORDENAMIENTO PARA IRFs (pi -> g -> h)
# El modelo principal es h -> g -> pi en todos los ejercicios de robustez.
# Para las graficas de IRF se usa el reordenamiento pi -> g -> h,
# identico al truco aplicado en pvar_estimacion.R, para que el efecto
# contemporaneo de h sobre pi sea visible en t=0.
# ============================================================
cat("\n============================================================\n")
cat("EJERCICIO 3: IRF CON REORDENAMIENTO (pi -> g -> h)\n")
cat("============================================================\n\n")

modelo_chol_alt <- estimar_pvar(panel, endogenas_irf, REZAGO_OPTIMO)

if (!is.null(modelo_chol_alt)) {
  irf_chol_alt <- calcular_irf(modelo_chol_alt, HORIZONTE_IRF, SEMILLA)
  graficar_irf_robustez(
    irf_chol_alt, endogenas_irf, "h_r", "pi_r",
    titulo      = "IRF Cholesky alternativo: inflacion ante choque de violencia",
    subtitulo   = paste0("Reordenamiento pi -> g -> h para IRF. Modelo: h -> g -> pi, p = ",
                         REZAGO_OPTIMO),
    ruta_salida = file.path(ruta_graficas, "g17_irf_cholesky_alt.png")
  )
}

# ============================================================
# EJERCICIO 4: CONTROL DE MERCANCIAS ALIMENTICIAS NACIONALES
# Se agrega el subindice de mercancias alimenticias del INPC
# (SF43405, Banxico) como exogena adicional en el partialling out.
# El objetivo es absorber la variacion nacional de precios de
# alimentos y verificar si el nulo persiste en la variacion
# regional residual, que es donde operaria el canal logistico
# de extorsion y cobro de piso.
# ============================================================
cat("\n============================================================\n")
cat("EJERCICIO 4: CONTROL DE MERCANCIAS ALIMENTICIAS NACIONALES\n")
cat("============================================================\n\n")

setToken(BANXICO_TOKEN)

resp_alim <- getSeriesData("SF43405",
                           startDate = FECHA_INICIO,
                           endDate   = FECHA_FIN)

alim_nacional <- as.data.frame(resp_alim[["SF43405"]]) %>%
  rename(fecha = date, inpc_alim_nac = value) %>%
  mutate(
    fecha      = as_date(fecha),
    d_alim_nac = log(inpc_alim_nac) - lag(log(inpc_alim_nac))
  ) %>%
  dplyr::select(fecha, d_alim_nac) %>%
  filter(!is.na(d_alim_nac))

panel_alim <- panel %>%
  left_join(alim_nacional, by = "fecha")

exogenas_alim <- c(exogenas_base, "d_alim_nac")

modelo_alim <- estimar_pvar(panel_alim, endogenas_principal,
                            REZAGO_OPTIMO, exogenas = exogenas_alim)

if (!is.null(modelo_alim)) {
  cat("Causalidad de Granger con control de alimentos:\n")
  granger_alim <- rbind(
    probar_granger(modelo_alim, "h_r", "pi_r"),
    probar_granger(modelo_alim, "h_r", "g_r"),
    probar_granger(modelo_alim, "g_r", "pi_r")
  )
  print(granger_alim)
  
  modelo_alim_irf <- estimar_pvar(panel_alim, endogenas_irf,
                                  REZAGO_OPTIMO, exogenas = exogenas_alim)
  
  if (!is.null(modelo_alim_irf)) {
    irf_alim <- calcular_irf(modelo_alim_irf, HORIZONTE_IRF, SEMILLA)
    graficar_irf_robustez(
      irf_alim, endogenas_irf, "h_r", "pi_r",
      titulo      = "IRF con control de mercancias alimenticias: inflacion ante choque de violencia",
      subtitulo   = paste0("Exogena adicional: d_alim_nac (SF43405 Banxico). p = ", REZAGO_OPTIMO),
      ruta_salida = file.path(ruta_graficas, "g18_irf_control_alimentos.png")
    )
  }
}

# ============================================================
# EJERCICIO 5: SEGUNDA ETAPA - VAR POR ESTADO
# Para cada estado se estima un VAR(p) individual con h, g, pi.
# Se extrae el IRF acumulado de pi ante un choque en h a
# HORIZONTE_ETAPA2 meses como escalar de sensibilidad estatal.
# Ese escalar se regresa en segunda etapa sobre caracteristicas
# estructurales del estado que proxiean la vulnerabilidad
# logistica al crimen organizado.
#
# PENDIENTE: completar el data frame "caracteristicas" con datos
# reales del Anuario Estadistico y Geografico de INEGI 2023
# antes de interpretar la regresion de segunda etapa.
# ============================================================
cat("\n============================================================\n")
cat("EJERCICIO 5: SEGUNDA ETAPA - HETEROGENEIDAD ESTRUCTURAL\n")
cat("============================================================\n\n")

# Caracteristicas estructurales por estado.
# pib_primario: participacion del sector primario en el PIB estatal.
#   Proxy de dependencia de cadenas logisticas expuestas a extorsion.
# densidad_vial: km de red carretera por km2 de superficie estatal.
#   Proxy inverso de aislamiento logistico (menor densidad = mayor vulnerabilidad).
caracteristicas <- data.frame(
  cve_ent       = 1:32,
  pib_primario  = rep(NA_real_, 32),
  densidad_vial = rep(NA_real_, 32)
)

cat("PENDIENTE: completar el data frame 'caracteristicas' con datos reales de INEGI.\n\n")

# Estima un VAR individual para un estado y devuelve el IRF acumulado
# de pi ante un choque en h a lo largo de "horizonte" meses.
# Devuelve NA si el panel del estado es demasiado corto para estimar.
extraer_irf_acumulado <- function(datos_estado, rezagos, horizonte) {
  vars_modelo <- datos_estado %>%
    dplyr::select(h, g, pi) %>%
    filter(complete.cases(.))
  
  if (nrow(vars_modelo) < rezagos * 3 + 20) return(NA_real_)
  
  modelo_var <- tryCatch(
    VAR(vars_modelo, p = rezagos, type = "const"),
    error = function(e) NULL
  )
  if (is.null(modelo_var)) return(NA_real_)
  
  irf_var <- tryCatch(
    irf(modelo_var, impulse = "h", response = "pi",
        n.ahead = horizonte, ortho = TRUE, boot = FALSE),
    error = function(e) NULL
  )
  if (is.null(irf_var)) return(NA_real_)
  
  sum(irf_var$irf$h[, "pi"], na.rm = TRUE)
}

estados_unicos <- sort(unique(panel$cve_ent))

irf_por_estado <- sapply(estados_unicos, function(estado) {
  datos_e <- panel %>%
    filter(cve_ent == estado) %>%
    arrange(periodo)
  extraer_irf_acumulado(datos_e, REZAGO_OPTIMO, HORIZONTE_ETAPA2)
})

resultados_segunda_etapa <- data.frame(
  cve_ent  = estados_unicos,
  irf_acum = irf_por_estado
) %>%
  left_join(caracteristicas, by = "cve_ent") %>%
  filter(!is.na(irf_acum))

cat("IRF acumulado h -> pi por estado (", HORIZONTE_ETAPA2, "meses):\n")
print(resultados_segunda_etapa %>%
        dplyr::select(cve_ent, irf_acum) %>%
        arrange(desc(irf_acum)), n = 32)

if (!all(is.na(caracteristicas$pib_primario))) {
  reg_segunda <- lm(irf_acum ~ pib_primario + densidad_vial,
                    data = resultados_segunda_etapa)
  cat("\nRegresion de segunda etapa:\n")
  print(summary(reg_segunda))
} else {
  cat("\nRegresion de segunda etapa pendiente: completar datos estructurales.\n")
}

write_csv(
  resultados_segunda_etapa,
  file.path(RUTA_BASE, "output", "irf_segunda_etapa.csv")
)

cat("IRF por estado exportado a: irf_segunda_etapa.csv\n")

# ============================================================
# EXPORTACION DEL RESUMEN
# ============================================================
sink(ruta_resumen)

cat("============================================================\n")
cat("RESUMEN DE ROBUSTEZ\n")
cat("Panel: 32 estados, enero 2004 - diciembre 2024\n")
cat("Rezago utilizado:", REZAGO_OPTIMO, "(MMSC-BIC de pvar_estimacion.R)\n")
cat("Umbral de violencia:", UMBRAL_VIOLENCIA, "por 100k (tasa estatal-mes promedio)\n")
cat("============================================================\n\n")

cat("Clasificacion de estados:\n")
print(tasa_promedio_estado %>%
        mutate(grupo = ifelse(tasa_media >= UMBRAL_VIOLENCIA,
                              "violento", "pacifico")) %>%
        arrange(desc(tasa_media)), n = 32)

cat("\nEstados violentos:", length(estados_violentos), "\n")
cat("Estados pacificos:", length(estados_pacificos), "\n\n")

cat("------------------------------------------------------------\n")
cat("EJERCICIO 1: ESTADOS VIOLENTOS\n")
cat("------------------------------------------------------------\n")
if (!is.null(modelo_violentos)) {
  print(summary(modelo_violentos))
} else {
  cat("Modelo no estimado (subpanel insuficiente).\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 1: ESTADOS PACIFICOS\n")
cat("------------------------------------------------------------\n")
if (!is.null(modelo_pacificos)) {
  print(summary(modelo_pacificos))
} else {
  cat("Modelo no estimado.\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 2: SUB-PERIODO 2004-2012\n")
cat("------------------------------------------------------------\n")
cat("Observaciones:", nrow(panel_sp1), "\n")
if (!is.null(modelo_sp1)) {
  print(summary(modelo_sp1))
} else {
  cat("Modelo no estimado.\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 2: SUB-PERIODO 2013-2024\n")
cat("------------------------------------------------------------\n")
cat("Observaciones:", nrow(panel_sp2), "\n")
if (!is.null(modelo_sp2)) {
  print(summary(modelo_sp2))
} else {
  cat("Modelo no estimado.\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 3: IRF CON REORDENAMIENTO (pi -> g -> h)\n")
cat("------------------------------------------------------------\n")
if (!is.null(modelo_chol_alt)) {
  print(summary(modelo_chol_alt))
} else {
  cat("Modelo no estimado.\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 4: CONTROL DE MERCANCIAS ALIMENTICIAS\n")
cat("------------------------------------------------------------\n")
if (!is.null(modelo_alim)) {
  print(summary(modelo_alim))
  cat("\nCausalidad de Granger:\n")
  print(granger_alim)
} else {
  cat("Modelo no estimado.\n")
}

cat("\n------------------------------------------------------------\n")
cat("EJERCICIO 5: IRF ACUMULADO POR ESTADO (segunda etapa)\n")
cat("------------------------------------------------------------\n")
cat("Horizonte acumulado:", HORIZONTE_ETAPA2, "meses\n")
print(resultados_segunda_etapa %>%
        dplyr::select(cve_ent, irf_acum) %>%
        arrange(desc(irf_acum)), n = 32)

sink()

cat("\nResumen exportado a:", ruta_resumen, "\n")