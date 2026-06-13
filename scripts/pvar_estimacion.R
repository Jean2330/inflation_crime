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
library(tidyr)

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_panel    <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")
ruta_resumen  <- file.path(RUTA_BASE, "output", "pvar_resultados.txt")
ruta_graficas <- file.path(RUTA_BASE, "output", "graficas")

dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# PARAMETROS
# ============================================================
REZAGOS_MAX     <- 3
HORIZONTE_IRF   <- 24
N_BOOT          <- 200   # replicas para bootstrap por bloques
SEMILLA         <- 2024
NIVEL_CONFIANZA <- 0.90  # bandas al 90%

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

cat("Panel:", nrow(panel), "obs |",
    n_distinct(panel$cve_ent), "estados |",
    n_distinct(panel$periodo), "periodos\n\n")

# ============================================================
# PARTIALLING OUT DE LAS EXOGENAS
# ============================================================
nombres_dummies <- paste0("d_mes", str_pad(2:12, 2, pad = "0"))
exogenas        <- c("d_tc", "cetes28", "expectativa", "d_wti", nombres_dummies)

cat("Aplicando partialling out de exogenas...\n")

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

# Aplica partialling out a un panel dado. Devuelve el panel con h_r, g_r, pi_r.
aplicar_partialling <- function(datos, vars_exog = exogenas) {
  datos %>%
    arrange(cve_ent, periodo) %>%
    group_by(cve_ent) %>%
    group_modify(function(d, key) {
      d$h_r  <- proyectar_serie(d$h,  d, vars_exog)
      d$g_r  <- proyectar_serie(d$g,  d, vars_exog)
      d$pi_r <- proyectar_serie(d$pi, d, vars_exog)
      d
    }) %>%
    ungroup() %>%
    as.data.frame()
}

panel_limpio <- aplicar_partialling(panel)

cat("NAs en h_r:", sum(is.na(panel_limpio$h_r)),
    "| g_r:", sum(is.na(panel_limpio$g_r)),
    "| pi_r:", sum(is.na(panel_limpio$pi_r)), "\n\n")

# ============================================================
# SELECCION DE REZAGOS POR MMSC
# ============================================================
cat("============================================================\n")
cat("SELECCION DE REZAGOS POR MMSC\n")
cat("============================================================\n\n")

set.seed(SEMILLA)
tabla_mmsc <- data.frame(rezagos  = integer(0),
                         mmsc_aic = numeric(0),
                         mmsc_bic = numeric(0),
                         mmsc_hqc = numeric(0))

for (p in 1:REZAGOS_MAX) {
  cat("Estimando p =", p, "...\n")
  modelo_p <- tryCatch(
    pvargmm(
      dependent_vars   = c("h_r", "g_r", "pi_r"),
      lags             = p,
      transformation   = "fod",
      data             = panel_limpio,
      panel_identifier = c(1, 2),
      steps            = "onestep",
      collapse         = TRUE
    ),
    error = function(e) {
      cat("  Error en p =", p, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(modelo_p)) {
    crit <- tryCatch(Andrews_Lu_MMSC(modelo_p),
                     error = function(e) {
                       cat("  Andrews_Lu_MMSC error:", conditionMessage(e), "\n")
                       NULL
                     })
    if (!is.null(crit)) {
      aic_val <- as.numeric(crit[["MMSC_AIC"]][[1]])
      bic_val <- as.numeric(crit[["MMSC_BIC"]][[1]])
      hqc_val <- as.numeric(crit[["MMSC_HQIC"]][[1]])
      cat(sprintf("  MMSC-AIC: %.2f | MMSC-BIC: %.2f | MMSC-HQC: %.2f\n",
                  aic_val, bic_val, hqc_val))
      tabla_mmsc <- rbind(tabla_mmsc,
                          data.frame(rezagos  = p,
                                     mmsc_aic = aic_val,
                                     mmsc_bic = bic_val,
                                     mmsc_hqc = hqc_val))
    } else {
      tabla_mmsc <- rbind(tabla_mmsc,
                          data.frame(rezagos  = p,
                                     mmsc_aic = NA_real_,
                                     mmsc_bic = NA_real_,
                                     mmsc_hqc = NA_real_))
    }
  }
}

cat("\nTabla MMSC:\n")
print(tabla_mmsc)

if (nrow(tabla_mmsc) == 0) stop("Ningun rezago pudo estimarse.")

if (all(is.na(tabla_mmsc$mmsc_bic))) {
  cat("MMSC-BIC no disponible. Seleccionando p = 1 por parsimonia.\n")
  rezagos_optimos <- 1L
} else {
  rezagos_optimos <- tabla_mmsc$rezagos[which.max(tabla_mmsc$mmsc_bic)]
}
cat("\nRezago optimo seleccionado:", rezagos_optimos, "\n\n")

# ============================================================
# ESTIMACION DEL MODELO PRINCIPAL
# ============================================================
cat("============================================================\n")
cat("ESTIMACION PVAR-GMM (p =", rezagos_optimos, ")\n")
cat("============================================================\n\n")

set.seed(SEMILLA)

modelo_pvar <- pvargmm(
  dependent_vars   = c("h_r", "g_r", "pi_r"),
  lags             = rezagos_optimos,
  transformation   = "fod",
  data             = panel_limpio,
  panel_identifier = c(1, 2),
  steps            = "onestep",
  collapse         = TRUE
)

print(summary(modelo_pvar))

# ============================================================
# VERIFICACION DE ESTABILIDAD
# ============================================================
cat("============================================================\n")
cat("VERIFICACION DE ESTABILIDAD\n")
cat("============================================================\n\n")

estabilidad <- tryCatch(stability(modelo_pvar),
                        error = function(e)
                          tryCatch(get_var_stability(modelo_pvar),
                                   error = function(e2) NULL))

if (!is.null(estabilidad)) {
  modulos <- estabilidad$Modulus
  cat("Modulos de los valores propios (descendente):\n")
  print(round(sort(modulos, decreasing = TRUE), 6))
  if (all(modulos < 1)) {
    cat("El modelo es ESTABLE (todos los modulos < 1).\n\n")
  } else {
    cat("ADVERTENCIA: el modelo puede ser INESTABLE.\n\n")
  }
} else {
  cat("No se pudo calcular estabilidad en esta version del paquete.\n\n")
  modulos <- numeric(0)
}

# ============================================================
# CAUSALIDAD DE GRANGER EN PANEL
# ============================================================
probar_granger <- function(modelo, causa, causada) {
  vars        <- modelo$dependent_vars
  n_vars      <- length(vars)
  n_lags      <- modelo$lags
  col_causa   <- which(vars == causa)
  col_causada <- which(vars == causada)
  if (length(col_causa) == 0 || length(col_causada) == 0) return(NULL)
  
  idx_coef <- sapply(1:n_lags, function(l) (l - 1) * n_vars + col_causa)
  beta_sub <- coef(modelo)[col_causada, idx_coef]
  se_sub   <- modelo$standard_error_first_step[col_causada, idx_coef]
  W_stat   <- tryCatch(sum((beta_sub / se_sub)^2), error = function(e) NA_real_)
  p_val    <- 1 - pchisq(W_stat, df = n_lags)
  cat(sprintf("  %s -> %s: Wald = %.3f, df = %d, p = %.4f\n",
              causa, causada, W_stat, n_lags, p_val))
  data.frame(causa = causa, causada = causada,
             wald = W_stat, df = n_lags, p_valor = p_val)
}

cat("Pruebas de causalidad de Granger (Wald):\n\n")
granger_resultado <- rbind(
  probar_granger(modelo_pvar, "h_r",  "pi_r"),
  probar_granger(modelo_pvar, "h_r",  "g_r"),
  probar_granger(modelo_pvar, "g_r",  "pi_r"),
  probar_granger(modelo_pvar, "pi_r", "h_r")
)
print(granger_resultado)

# ============================================================
# BOOTSTRAP POR BLOQUES PARA BANDAS DE CONFIANZA DE LAS IRF
#
# Se remuestrean los N grupos (estados) con reemplazo. Para cada
# muestra bootstrap se re-aplica el partialling out, se reestima
# el PVAR con los mismos hiperparametros, y se calcula la IRF
# ortogonalizada. Las bandas son los percentiles (alfa/2) y
# (1 - alfa/2) de la distribucion empirica de las IRF.
#
# Este procedimiento es equivalente al bootstrap de panel por
# grupos descrito en Gonçalves (2011) y es computacionalmente
# viable porque cada replica usa N = 32 grupos, no el panel
# completo sin remuestreo.
# ============================================================
cat("\n============================================================\n")
cat("BOOTSTRAP POR BLOQUES PARA BANDAS IRF (N_BOOT =", N_BOOT, ")\n")
cat("============================================================\n\n")

estados_unicos <- sort(unique(panel_limpio$cve_ent))

# Calcula la IRF puntual de un modelo estimado.
# Devuelve un vector de longitud (horizonte + 1) con la respuesta
# de la variable "respuesta" ante un choque en "impulso".
calcular_irf_puntual <- function(modelo, endogenas, impulso, respuesta, horizonte) {
  col_impulso <- which(endogenas == impulso)
  irf_obj     <- tryCatch(oirf(modelo, n.ahead = horizonte),
                          error = function(e) NULL)
  if (is.null(irf_obj) || !respuesta %in% names(irf_obj)) return(NULL)
  irf_obj[[respuesta]][, col_impulso]
}

# Una replica del bootstrap: remuestrea grupos, reestima, devuelve IRF.
una_replica_boot <- function(panel_base, endogenas, rezagos, horizonte,
                             impulso, respuesta, vars_exog) {
  # Remuestrear estados con reemplazo
  estados_boot  <- sample(estados_unicos, size = length(estados_unicos),
                          replace = TRUE)
  # Construir panel bootstrap reasignando cve_ent para que sean unicos
  panel_boot <- lapply(seq_along(estados_boot), function(i) {
    panel_base %>%
      filter(cve_ent == estados_boot[i]) %>%
      mutate(cve_ent = i)
  })
  panel_boot <- do.call(rbind, panel_boot) %>%
    arrange(cve_ent, periodo) %>%
    as.data.frame()
  
  # Partialling out sobre la muestra bootstrap
  panel_boot_limpio <- tryCatch(
    aplicar_partialling(panel_boot, vars_exog),
    error = function(e) NULL
  )
  if (is.null(panel_boot_limpio)) return(NULL)
  
  # Reestimar PVAR
  modelo_boot <- tryCatch(
    pvargmm(
      dependent_vars   = endogenas,
      lags             = rezagos,
      transformation   = "fod",
      data             = panel_boot_limpio,
      panel_identifier = c(1, 2),
      steps            = "onestep",
      collapse         = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(modelo_boot)) return(NULL)
  
  calcular_irf_puntual(modelo_boot, endogenas, impulso, respuesta, horizonte)
}

# Ejecuta el bootstrap y devuelve una lista con media, lower y upper
# para graficar directamente.
bootstrap_irf <- function(panel_base, endogenas, rezagos, horizonte,
                          impulso, respuesta, vars_exog,
                          n_boot, semilla, nivel) {
  set.seed(semilla)
  alfa <- 1 - nivel
  
  # IRF puntual sobre el modelo completo
  modelo_completo <- tryCatch(
    pvargmm(
      dependent_vars   = endogenas,
      lags             = rezagos,
      transformation   = "fod",
      data             = aplicar_partialling(panel_base, vars_exog),
      panel_identifier = c(1, 2),
      steps            = "onestep",
      collapse         = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(modelo_completo)) {
    cat("No se pudo estimar el modelo completo para bootstrap.\n")
    return(NULL)
  }
  
  irf_puntual <- calcular_irf_puntual(modelo_completo, endogenas,
                                      impulso, respuesta, horizonte)
  if (is.null(irf_puntual)) {
    cat("No se pudo calcular la IRF puntual.\n")
    return(NULL)
  }
  
  n_horizonte <- length(irf_puntual)
  
  # Acumular replicas
  replicas <- matrix(NA_real_, nrow = n_boot, ncol = n_horizonte)
  cat("Corriendo bootstrap")
  for (b in 1:n_boot) {
    if (b %% 50 == 0) cat(".", b)
    rep_b <- una_replica_boot(panel_base, endogenas, rezagos, horizonte,
                              impulso, respuesta, vars_exog)
    if (!is.null(rep_b) && length(rep_b) == n_horizonte) {
      replicas[b, ] <- rep_b
    }
  }
  cat("\n")
  
  # Descartar replicas fallidas
  replicas_ok <- replicas[complete.cases(replicas), ]
  cat("Replicas exitosas:", nrow(replicas_ok), "de", n_boot, "\n")
  
  if (nrow(replicas_ok) < 10) {
    cat("Replicas insuficientes. Se grafica sin bandas.\n")
    return(data.frame(
      horizonte = 0:(n_horizonte - 1),
      media     = irf_puntual,
      lower     = irf_puntual,
      upper     = irf_puntual
    ))
  }
  
  data.frame(
    horizonte = 0:(n_horizonte - 1),
    media     = irf_puntual,
    lower     = apply(replicas_ok, 2, quantile, probs = alfa / 2,
                      na.rm = TRUE),
    upper     = apply(replicas_ok, 2, quantile, probs = 1 - alfa / 2,
                      na.rm = TRUE)
  )
}

# Grafica una IRF con bandas de confianza.
graficar_irf <- function(datos_irf, titulo, subtitulo, ruta_salida) {
  if (is.null(datos_irf)) return(invisible(NULL))
  
  g <- ggplot(datos_irf, aes(x = horizonte)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = COLOR_PRINCIPAL, alpha = 0.2) +
    geom_line(aes(y = media), color = COLOR_PRINCIPAL, linewidth = 0.9) +
    geom_hline(yintercept = 0, color = COLOR_NEUTRO, linewidth = 0.5,
               linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, HORIZONTE_IRF, by = 6)) +
    labs(title    = titulo,
         subtitle = subtitulo,
         x = "Meses despues del choque", y = "Respuesta") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(ruta_salida, g, width = 8, height = 5, dpi = 150)
  cat("Grafica guardada:", ruta_salida, "\n")
}

# Ordenamiento pi -> g -> h para recuperar impacto contemporaneo en t = 0.
endogenas_irf <- c("pi_r", "g_r", "h_r")

subtitulo_base <- paste0("PVAR-GMM, p = ", rezagos_optimos,
                         ". Cholesky: h -> g -> pi. IRF ortogonalizada. ",
                         "Bandas al ", round(NIVEL_CONFIANZA * 100), "%")

# IRF principal: respuesta de pi ante choque en h
cat("\nCalculando IRF con bootstrap: pi ante h...\n")
datos_irf_pi_h <- bootstrap_irf(
  panel_base  = panel,
  endogenas   = endogenas_irf,
  rezagos     = rezagos_optimos,
  horizonte   = HORIZONTE_IRF,
  impulso     = "h_r",
  respuesta   = "pi_r",
  vars_exog   = exogenas,
  n_boot      = N_BOOT,
  semilla     = SEMILLA,
  nivel       = NIVEL_CONFIANZA
)

graficar_irf(
  datos_irf_pi_h,
  titulo      = "Respuesta de la inflacion ante un choque de violencia",
  subtitulo   = subtitulo_base,
  ruta_salida = file.path(ruta_graficas, "g10_irf_h_a_pi.png")
)

# IRF: respuesta de g ante choque en h
cat("\nCalculando IRF con bootstrap: g ante h...\n")
datos_irf_g_h <- bootstrap_irf(
  panel_base  = panel,
  endogenas   = endogenas_irf,
  rezagos     = rezagos_optimos,
  horizonte   = HORIZONTE_IRF,
  impulso     = "h_r",
  respuesta   = "g_r",
  vars_exog   = exogenas,
  n_boot      = N_BOOT,
  semilla     = SEMILLA,
  nivel       = NIVEL_CONFIANZA
)

graficar_irf(
  datos_irf_g_h,
  titulo      = "Respuesta de la brecha del producto ante un choque de violencia",
  subtitulo   = subtitulo_base,
  ruta_salida = file.path(ruta_graficas, "g11_irf_h_a_g.png")
)

# ============================================================
# FEVD
# ============================================================
fevd_resultado <- fevd_orthogonal(modelo_pvar, n.ahead = HORIZONTE_IRF)

cat("FEVD de pi_r (proporcion explicada por cada variable):\n")
print(round(fevd_resultado$pi_r, 4))

fevd_pi <- as.data.frame(fevd_resultado$pi_r)
fevd_pi$horizonte <- 1:nrow(fevd_pi)

fevd_largo <- pivot_longer(fevd_pi, cols = -horizonte,
                           names_to = "variable", values_to = "proporcion")

paleta_fevd    <- c(h_r = COLOR_ACENTO, g_r = COLOR_NEUTRO, pi_r = COLOR_PRINCIPAL)
etiquetas_fevd <- c(h_r = "Violencia (h)", g_r = "Brecha del producto (g)",
                    pi_r = "Inflacion (pi)")

g_fevd <- ggplot(fevd_largo,
                 aes(x = horizonte, y = proporcion, fill = variable)) +
  geom_area(alpha = 0.85, position = "stack") +
  scale_fill_manual(values = paleta_fevd, labels = etiquetas_fevd) +
  scale_x_continuous(breaks = seq(0, HORIZONTE_IRF, by = 6)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title    = "Descomposicion de varianza del error de prediccion de la inflacion",
       subtitle = paste0("PVAR-GMM, p = ", rezagos_optimos,
                         ". Horizonte ", HORIZONTE_IRF, " meses"),
       x = "Horizonte (meses)", y = "Proporcion de la varianza", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(ruta_graficas, "g12_fevd.png"),
       g_fevd, width = 9, height = 6, dpi = 150)
cat("FEVD guardado en g12_fevd.png\n")

# ============================================================
# EXPORTACION DEL RESUMEN
# ============================================================
sink(ruta_resumen)

cat("============================================================\n")
cat("RESULTADOS PVAR-GMM\n")
cat("Panel: 32 estados, enero 2004 - diciembre 2024\n")
cat("Vector endogeno: h -> g -> pi (Cholesky)\n")
cat("Exogenas proyectadas por partialling out:\n")
cat("  d_tc, cetes28, expectativa, d_wti, dummies de mes\n")
cat("Estimacion: GMM un paso, collapse = TRUE, fod (forward orthogonal deviations)\n")
cat("Bootstrap por bloques:", N_BOOT, "replicas. Nivel de confianza:",
    round(NIVEL_CONFIANZA * 100), "%\n")
cat("============================================================\n\n")

cat("Tabla MMSC:\n")
print(tabla_mmsc)
cat("\nRezago optimo (MMSC-BIC):", rezagos_optimos, "\n\n")

cat("Resumen del modelo:\n")
print(summary(modelo_pvar))

cat("\nEstabilidad (modulos de valores propios):\n")
if (length(modulos) > 0) {
  print(round(sort(modulos, decreasing = TRUE), 6))
  cat("Sistema estable:",
      ifelse(all(modulos < 1), "SI (todos < 1)", "NO"), "\n")
} else {
  cat("No disponible.\n")
}

cat("\nCausalidad de Granger (Wald):\n")
print(granger_resultado)

cat("\nFEVD de pi (horizonte", HORIZONTE_IRF, "meses):\n")
print(round(fevd_resultado$pi_r, 4))

sink()

cat("\nResumen exportado a:", ruta_resumen, "\n")