rm(list = ls())
cat("\014")

# ============================================================
#                       LIBRERIAS
# ============================================================
library(tidyverse)
library(lubridate)
library(scales)

# ============================================================
#                    CREDENCIALES Y PARAMS
# ============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))

RUTA_BASE <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

ruta_inpc    <- file.path(RUTA_BASE, "data", "raw", "inpc_ciudades.csv")
ruta_itaee   <- file.path(RUTA_BASE, "data", "raw", "itaee.csv")
ruta_narco   <- file.path(RUTA_BASE, "data", "raw", "narco_homicidios.csv")
ruta_macro   <- file.path(RUTA_BASE, "data", "raw", "macro_controles.csv")
ruta_resumen <- file.path(RUTA_BASE, "output", "diagnosticos_resumen.txt")
ruta_graficas <- file.path(RUTA_BASE, "output", "graficas")

dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

COLOR_PRINCIPAL <- "#2C5F8A"
COLOR_ACENTO    <- "#C94040"
COLOR_NEUTRO    <- "#A8A8A8"

# ============================================================
# BLOQUE 1: INPC regional
# ============================================================
cat("Cargando INPC...\n")

inpc <- read_csv(ruta_inpc, show_col_types = FALSE,
                 col_types = cols(cve_ent = col_integer())) %>%
  mutate(fecha = as_date(fecha)) %>%
  arrange(cve_ent, fecha)

cat("INPC:", nrow(inpc), "obs | estados:", n_distinct(inpc$cve_ent),
    "| meses:", n_distinct(inpc$fecha), "\n")
cat("Fechas:", format(min(inpc$fecha)), "a", format(max(inpc$fecha)), "\n")
cat("NAs en inpc:", sum(is.na(inpc$inpc)), "\n\n")

inpc_pi <- inpc %>%
  group_by(cve_ent, entidad) %>%
  mutate(pi = log(inpc) - lag(log(inpc))) %>%
  ungroup() %>%
  filter(!is.na(pi))

# Grafica 1: Inflacion promedio anual por estado
inpc_resumen_estado <- inpc_pi %>%
  group_by(entidad) %>%
  summarise(pi_anual = mean(pi, na.rm = TRUE) * 12 * 100)

g1 <- ggplot(inpc_resumen_estado,
             aes(x = pi_anual, y = reorder(entidad, pi_anual))) +
  geom_col(fill = COLOR_PRINCIPAL, width = 0.7) +
  geom_vline(xintercept = mean(inpc_resumen_estado$pi_anual),
             linetype = "dashed", color = COLOR_ACENTO, linewidth = 0.7) +
  labs(title    = "Inflacion promedio anual por estado, 2004-2024",
       subtitle = "Linea punteada = promedio nacional. Delta log mensual anualizado",
       x = "Inflacion anual promedio (%)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(ruta_graficas, "g01_inflacion_por_estado.png"),
       g1, width = 9, height = 8, dpi = 150)

# Grafica 2: Series de tiempo, estados seleccionados
inpc_pi_seleccion <- inpc_pi %>%
  filter(cve_ent %in% c(2, 5, 8, 9, 14, 19, 25, 26, 28, 31))

g2 <- ggplot(inpc_pi_seleccion, aes(x = fecha, y = pi * 100)) +
  geom_line(color = COLOR_PRINCIPAL, linewidth = 0.4) +
  geom_hline(yintercept = 0, color = COLOR_NEUTRO, linewidth = 0.3) +
  facet_wrap(~ entidad, ncol = 2, scales = "free_y") +
  labs(title    = "Inflacion mensual en estados seleccionados, 2004-2024",
       subtitle = "Delta log del INPC por ciudad, en puntos porcentuales",
       x = NULL, y = "Inflacion mensual (%)") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(ruta_graficas, "g02_inflacion_series_estados.png"),
       g2, width = 10, height = 12, dpi = 150)

# Grafica 3: Mapa de calor del panel de inflacion
inpc_heat <- inpc_pi %>%
  mutate(pi_clipped = pmin(pmax(pi * 100, -1.5), 1.5))

g3 <- ggplot(inpc_heat,
             aes(x = fecha, y = reorder(entidad, cve_ent), fill = pi_clipped)) +
  geom_tile() +
  scale_fill_gradient2(low = COLOR_ACENTO, mid = "white", high = COLOR_PRINCIPAL,
                       midpoint = 0, name = "pi (%)", limits = c(-1.5, 1.5)) +
  labs(title    = "Panel de inflacion mensual: mapa de calor",
       subtitle = "Valores recortados a [-1.5, 1.5] pp. Rojo = deflacion, azul = inflacion",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")

ggsave(file.path(ruta_graficas, "g03_heatmap_inflacion.png"),
       g3, width = 12, height = 8, dpi = 150)

# ============================================================
# BLOQUE 2: Narco-homicidios
# ============================================================
cat("Cargando narco-homicidios...\n")

narco <- read_csv(ruta_narco, show_col_types = FALSE) %>%
  mutate(
    cve_ent = as.integer(str_sub(str_pad(cve_inegi, 5, pad = "0"), 1, 2)),
    fecha   = dmy(fecha)
  ) %>%
  group_by(cve_ent, fecha) %>%
  summarise(
    homicidios = sum(homicidios_narco, na.rm = TRUE),
    poblacion  = sum(poblacion_total,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(tasa = homicidios / poblacion * 100000,
         h    = log(1 + tasa))

cat("Narco:", nrow(narco), "obs estado-mes | estados:", n_distinct(narco$cve_ent), "\n")
cat("Total homicidios narco:", sum(narco$homicidios), "\n")
cat("Tasa promedio nacional (por 100k):", round(mean(narco$tasa), 3), "\n\n")

# Grafica 4: Tasa nacional en el tiempo
narco_nacional <- narco %>%
  group_by(fecha) %>%
  summarise(homicidios = sum(homicidios), poblacion = sum(poblacion),
            tasa = homicidios / poblacion * 100000)

g4 <- ggplot(narco_nacional, aes(x = fecha, y = tasa)) +
  geom_line(color = COLOR_PRINCIPAL, linewidth = 0.7) +
  geom_area(fill = COLOR_PRINCIPAL, alpha = 0.15) +
  labs(title    = "Tasa nacional de narco-homicidios, 2004-2024",
       subtitle = "Homicidios por organizacion criminal por cada 100,000 habitantes",
       x = NULL, y = "Tasa por 100,000 hab.") +
  theme_minimal(base_size = 11)

ggsave(file.path(ruta_graficas, "g04_narco_nacional.png"),
       g4, width = 10, height = 5, dpi = 150)

# Grafica 5: Promedio por estado
narco_estado <- narco %>%
  group_by(cve_ent) %>%
  summarise(tasa_media = mean(tasa, na.rm = TRUE)) %>%
  left_join(distinct(select(inpc, cve_ent, entidad)), by = "cve_ent")

g5 <- ggplot(narco_estado,
             aes(x = tasa_media, y = reorder(entidad, tasa_media))) +
  geom_col(fill = COLOR_ACENTO, width = 0.7) +
  labs(title    = "Tasa promedio de narco-homicidios por estado, 2004-2024",
       subtitle = "Promedio mensual de homicidios por 100,000 habitantes",
       x = "Tasa promedio (por 100,000 hab.)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(ruta_graficas, "g05_narco_por_estado.png"),
       g5, width = 9, height = 8, dpi = 150)

# Grafica 6: Mapa de calor narco
narco_heat <- narco %>%
  left_join(distinct(select(inpc, cve_ent, entidad)), by = "cve_ent") %>%
  mutate(tasa_clipped = pmin(tasa, 15))

g6 <- ggplot(narco_heat,
             aes(x = fecha, y = reorder(entidad, cve_ent), fill = tasa_clipped)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = COLOR_ACENTO,
                      name = "Tasa\n(por 100k)", limits = c(0, 15)) +
  labs(title    = "Panel de narco-homicidios: mapa de calor",
       subtitle = "Tasa recortada en 15 por 100,000 hab. para legibilidad",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")

ggsave(file.path(ruta_graficas, "g06_heatmap_narco.png"),
       g6, width = 12, height = 8, dpi = 150)

# ============================================================
# BLOQUE 3: ITAEE
# ============================================================
cat("Cargando ITAEE...\n")

itaee <- read_csv(ruta_itaee, show_col_types = FALSE,
                  col_types = cols(cve_ent = col_integer())) %>%
  mutate(fecha = as_date(fecha)) %>%
  arrange(cve_ent, fecha)

cat("ITAEE:", nrow(itaee), "obs | NAs:", sum(is.na(itaee$itaee)), "\n\n")

itaee_crecimiento <- itaee %>%
  group_by(cve_ent, entidad) %>%
  mutate(g_itaee = log(itaee) - lag(log(itaee), 12)) %>%
  ungroup() %>%
  filter(!is.na(g_itaee))

# Grafica 7: Crecimiento anual ITAEE, estados seleccionados
g7 <- ggplot(filter(itaee_crecimiento, cve_ent %in% c(2, 5, 8, 9, 14, 19, 25, 26, 28, 31)),
             aes(x = fecha, y = g_itaee * 100)) +
  geom_line(color = COLOR_PRINCIPAL, linewidth = 0.4) +
  geom_hline(yintercept = 0, color = COLOR_NEUTRO, linewidth = 0.3) +
  facet_wrap(~ entidad, ncol = 2, scales = "free_y") +
  labs(title    = "Crecimiento anual del ITAEE en estados seleccionados",
       subtitle = "Log-diferencia a 12 meses, en puntos porcentuales",
       x = NULL, y = "Crecimiento anual (%)") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(ruta_graficas, "g07_itaee_crecimiento.png"),
       g7, width = 10, height = 12, dpi = 150)

# ============================================================
# BLOQUE 4: Controles macroeconomicos
# ============================================================
cat("Cargando controles macro...\n")

macro <- read_csv(ruta_macro, show_col_types = FALSE) %>%
  mutate(fecha = as_date(fecha)) %>%
  arrange(fecha) %>%
  mutate(d_tc  = log(tc_fix)   - lag(log(tc_fix)),
         d_wti = log(petroleo) - lag(log(petroleo)))

cat("Macro:", nrow(macro), "obs | NAs por columna:\n")
print(colSums(is.na(macro)))
cat("\n")

# Grafica 8: Panel de controles macro
macro_long <- macro %>%
  select(fecha, d_tc, cetes28, expectativa, d_wti) %>%
  pivot_longer(-fecha, names_to = "variable", values_to = "valor") %>%
  mutate(variable = recode(variable,
                           d_tc        = "Depreciacion TC (log-dif)",
                           cetes28     = "CETES 28 dias (%)",
                           expectativa = "Expectativa inflacion anual (%)",
                           d_wti       = "Cambio precio WTI (log-dif)"
  ))

g8 <- ggplot(macro_long, aes(x = fecha, y = valor)) +
  geom_line(color = COLOR_PRINCIPAL, linewidth = 0.5) +
  geom_hline(yintercept = 0, color = COLOR_NEUTRO, linewidth = 0.3) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  labs(title = "Variables de control macroeconomico, 2004-2024",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 9))

ggsave(file.path(ruta_graficas, "g08_controles_macro.png"),
       g8, width = 10, height = 10, dpi = 150)

# ============================================================
# BLOQUE 5: Correlacion preliminar violencia e inflacion
# ============================================================
narco_media <- narco %>%
  group_by(cve_ent) %>%
  summarise(h_media = mean(h, na.rm = TRUE))

inpc_media <- inpc_pi %>%
  group_by(cve_ent) %>%
  summarise(pi_media = mean(pi, na.rm = TRUE) * 100)

correlacion_cross <- narco_media %>%
  inner_join(inpc_media, by = "cve_ent") %>%
  left_join(distinct(select(inpc, cve_ent, entidad)), by = "cve_ent")

g9 <- ggplot(correlacion_cross, aes(x = h_media, y = pi_media)) +
  geom_point(color = COLOR_PRINCIPAL, size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = COLOR_ACENTO,
              fill = COLOR_ACENTO, alpha = 0.1, linewidth = 0.8) +
  geom_text(aes(label = entidad), size = 2.5, vjust = -0.7,
            color = "gray40", check_overlap = TRUE) +
  labs(title    = "Correlacion bruta entre violencia e inflacion por estado",
       subtitle = "Promedios 2004-2024. h = log(1 + tasa narco), pi = delta log INPC mensual",
       x = "Exposicion criminal promedio h = log(1 + tasa)",
       y = "Inflacion mensual promedio pi (%)") +
  theme_minimal(base_size = 11)

ggsave(file.path(ruta_graficas, "g09_correlacion_bruta.png"),
       g9, width = 9, height = 7, dpi = 150)

r_pearson <- cor(correlacion_cross$h_media, correlacion_cross$pi_media,
                 use = "complete.obs")
cat("Correlacion bruta cross-sectional (h, pi):", round(r_pearson, 4), "\n\n")

# ============================================================
# BLOQUE 6: Resumen al archivo de texto
# ============================================================
sink(ruta_resumen)

cat("============================================================\n")
cat("DIAGNOSTICOS Y ESTADISTICAS DESCRIPTIVAS\n")
cat("Panel: 32 estados x 252 meses, enero 2004 - diciembre 2024\n")
cat("============================================================\n\n")

cat("------ INPC (inflacion mensual pi) ------\n")
inpc_pi %>% select(pi) %>% mutate(pi = pi * 100) %>% summary() %>% print()

cat("\n------ Narco-homicidios (tasa por 100k) ------\n")
narco %>% select(tasa) %>% summary() %>% print()

cat("\nEstados con tasa promedio > 5 por 100k:\n")
narco_estado %>% filter(tasa_media > 5) %>% arrange(desc(tasa_media)) %>% print()

cat("\n------ ITAEE (nivel, interpolado a mensual) ------\n")
itaee %>% select(itaee) %>% summary() %>% print()

cat("\n------ Controles macro ------\n")
macro %>% select(d_tc, cetes28, expectativa, d_wti) %>% summary() %>% print()

cat("\n------ Correlacion bruta cross-sectional ------\n")
cat("Pearson (h_media, pi_media):", round(r_pearson, 4), "\n")
cat("Un valor bajo no descarta el efecto causal: la correlacion bruta\n")
cat("no condiciona por efectos fijos ni controles macroeconomicos.\n")

sink()

cat("Resumen guardado en:  ", ruta_resumen, "\n")
cat("Graficas guardadas en:", ruta_graficas, "\n")
cat("Total graficas generadas: 9\n")