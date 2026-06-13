rm(list = ls())
cat("\014")

# ============================================================
# LIBRERIAS
# ============================================================
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)

# ============================================================
# RUTAS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))
RUTA_BASE  <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
ruta_panel <- file.path(RUTA_BASE, "data", "processed", "panel_completo.csv")
ruta_salida <- file.path(RUTA_BASE, "output", "graficas", "g_series_nacionales.png")

# ============================================================
# COLORES
# ============================================================
COLOR_H  <- "#2C5F8A"   # azul principal para violencia
COLOR_G  <- "#2C5F8A"   # mismo azul para brecha
COLOR_PI <- "#2C5F8A"   # mismo azul para inflacion
COLOR_CERO <- "#A8A8A8" # gris para linea de cero

# ============================================================
# LECTURA Y AGREGACION NACIONAL
# ============================================================
cat("Cargando panel...\n")

panel <- read_csv(ruta_panel,
                  col_types = cols(cve_ent = col_integer()),
                  show_col_types = FALSE) %>%
  mutate(fecha = as_date(fecha))

# Promedios nacionales por mes (promedio simple entre estados)
nacional <- panel %>%
  group_by(fecha) %>%
  summarise(
    h_nacional  = mean(h,  na.rm = TRUE),
    g_nacional  = mean(g,  na.rm = TRUE),
    pi_nacional = mean(pi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(fecha)

cat("Periodos:", nrow(nacional), "\n")
cat("Rango:", format(min(nacional$fecha), "%Y-%m"),
    "a", format(max(nacional$fecha), "%Y-%m"), "\n\n")

# ============================================================
# TEMA BASE COMPARTIDO
# ============================================================
tema_base <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92"),
    axis.title.x      = element_blank(),
    plot.title        = element_text(size = 11, face = "plain"),
    plot.margin       = margin(t = 4, r = 8, b = 4, l = 8)
  )

# ============================================================
# PANEL SUPERIOR: tasa de narco-homicidios
# ============================================================
p_h <- ggplot(nacional, aes(x = fecha, y = h_nacional)) +
  geom_area(fill = COLOR_H, alpha = 0.18) +
  geom_line(color = COLOR_H, linewidth = 0.7) +
  geom_hline(yintercept = 0, color = COLOR_CERO,
             linewidth = 0.4, linetype = "dashed") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               limits = range(nacional$fecha)) +
  labs(
    title = expression(italic(h)[it] ~ ": tasa de narco-homicidios (promedio nacional)"),
    y = "log(1 + tasa)"
  ) +
  tema_base

# ============================================================
# PANEL MEDIO: brecha del producto
# ============================================================
p_g <- ggplot(nacional, aes(x = fecha, y = g_nacional)) +
  geom_hline(yintercept = 0, color = COLOR_CERO,
             linewidth = 0.4, linetype = "dashed") +
  geom_area(aes(y = ifelse(g_nacional >= 0, g_nacional, 0)),
            fill = COLOR_G, alpha = 0.18) +
  geom_area(aes(y = ifelse(g_nacional < 0, g_nacional, 0)),
            fill = "#C94040", alpha = 0.18) +
  geom_line(color = COLOR_G, linewidth = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               limits = range(nacional$fecha)) +
  labs(
    title = expression(italic(g)[it] ~ ": brecha del producto (promedio nacional)"),
    y = "Desviacion"
  ) +
  tema_base

# ============================================================
# PANEL INFERIOR: inflacion mensual
# ============================================================
p_pi <- ggplot(nacional, aes(x = fecha, y = pi_nacional * 100)) +
  geom_hline(yintercept = 0, color = COLOR_CERO,
             linewidth = 0.4, linetype = "dashed") +
  geom_line(color = COLOR_PI, linewidth = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               limits = range(nacional$fecha)) +
  labs(
    title = expression(italic(pi)[it] ~ ": inflacion mensual (promedio nacional)"),
    y = "Puntos porcentuales",
    x = NULL
  ) +
  tema_base

# ============================================================
# FIGURA COMBINADA
# ============================================================
figura <- arrangeGrob(
  p_h, p_g, p_pi,
  ncol   = 1,
  heights = c(1, 1, 1)
)

ggsave(
  filename = ruta_salida,
  plot     = figura,
  width    = 10,
  height   = 9,
  dpi      = 150
)

cat("Figura guardada en:", ruta_salida, "\n")