rm(list = ls())
cat("\014")

# =============================================================
# LIBRERIAS
# =============================================================
library(readr)
library(tidyverse)
library(lubridate)
library(here)

# =============================================================
# CREDENCIALES 
# =============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))

# =============================================================
# LECTURA
# =============================================================
inpc_raw <- read_csv(URL_INPC_CRUDO,
                     col_names = FALSE,
                     show_col_types = FALSE,
                     locale = locale(encoding = "latin1"))

# Busqueda automatica y limpieza de encabezados
fila_titulos <- which(str_detect(inpc_raw[[1]], "T.tulo"))[1]
titulos <- as.character(inpc_raw[fila_titulos, ])

nombres_col <- case_when(
  str_detect(titulos, "Resumen, Principales .ndices") ~ "Nacional",
  str_detect(titulos, "Por ciudad") ~ str_extract(titulos, "(?<=Por ciudad, ).*?(?=, .ndice)"),
  is.na(titulos) ~ "basura",
  TRUE ~ "periodo"
)

nombres_col[1] <- "periodo"
nombres_col <- str_remove(nombres_col, "^\\d+\\.\\s*")
nombres_col <- make.unique(nombres_col)

names(inpc_raw) <- nombres_col

inpc <- inpc_raw %>%
  filter(str_detect(periodo, "^[a-z]{3}-\\d{2}")) %>%
  select(-starts_with("basura"))

# =============================================================
#                     TRANSFORMACION
# =============================================================
meses_espanol <- c("ene" = 1, "feb" = 2, "mar" = 3, "abr" = 4, "may" = 5, "jun" = 6,
                   "jul" = 7, "ago" = 8, "sep" = 9, "oct" = 10, "nov" = 11, "dic" = 12)

inpc_largo <- inpc %>%
  pivot_longer(cols = -periodo, names_to = "ciudad", values_to = "inpc") %>%
  mutate(
    inpc     = suppressWarnings(as.numeric(inpc)),
    mes_str  = str_sub(periodo, 1, 3),
    anio_str = as.integer(str_sub(periodo, 5, 6)),
    anio     = ifelse(anio_str > 50, 1900 + anio_str, 2000 + anio_str),
    mes      = meses_espanol[mes_str],
    fecha    = make_date(year = anio, month = mes, day = 1)
  ) %>%
  drop_na(fecha, inpc) %>%
  select(ciudad, fecha, inpc) %>%
  arrange(ciudad, fecha)

# Duplicar el AMCM para asignarlo tambien al Estado de Mexico.
# El Area Metropolitana cubre municipios conurbados de CDMX y Edomex.
inpc_largo <- inpc_largo %>%
  bind_rows(
    inpc_largo %>%
      filter(ciudad == "Area Metropolitana de la Cd. de México") %>%
      mutate(ciudad = "Area Metropolitana de la Cd. de México [Edomex]")
  )

# Asignacion de claves de estado.
# Baja California usa (?!S) para no capturar Baja California Sur.
inpc_con_claves <- inpc_largo %>%
  mutate(
    cve_ent = case_when(
      str_detect(ciudad, "Ags\\.") ~ "01",
      str_detect(ciudad, "B\\.C\\.(?!S)") ~ "02",
      str_detect(ciudad, "B\\.C\\.S\\.") ~ "03",
      str_detect(ciudad, "Camp\\.") ~ "04",
      str_detect(ciudad, "Coah\\.") ~ "05",
      str_detect(ciudad, "Col\\.") ~ "06",
      str_detect(ciudad, "Chis\\.") ~ "07",
      str_detect(ciudad, "Chih\\.") ~ "08",
      str_detect(ciudad, "\\[Edomex\\]") ~ "15",
      str_detect(ciudad, "Cd\\. de México|Ciudad de México") ~ "09",
      str_detect(ciudad, "Dgo\\.") ~ "10",
      str_detect(ciudad, "Gto\\.") ~ "11",
      str_detect(ciudad, "Gro\\.") ~ "12",
      str_detect(ciudad, "Hgo\\.") ~ "13",
      str_detect(ciudad, "Jal\\.") ~ "14",
      str_detect(ciudad, "Mich\\.") ~ "16",
      str_detect(ciudad, "Mor\\.") ~ "17",
      str_detect(ciudad, "Nay\\.") ~ "18",
      str_detect(ciudad, "N\\.L\\.") ~ "19",
      str_detect(ciudad, "Oax\\.") ~ "20",
      str_detect(ciudad, "Pue\\.") ~ "21",
      str_detect(ciudad, "Qro\\.") ~ "22",
      str_detect(ciudad, "Q\\. Roo|Q\\.R\\.") ~ "23",
      str_detect(ciudad, "S\\.L\\.P\\.") ~ "24",
      str_detect(ciudad, "Sin\\.") ~ "25",
      str_detect(ciudad, "Son\\.") ~ "26",
      str_detect(ciudad, "Tab\\.") ~ "27",
      str_detect(ciudad, "Tamps\\.") ~ "28",
      str_detect(ciudad, "Tlax\\.") ~ "29",
      str_detect(ciudad, "Ver\\.") ~ "30",
      str_detect(ciudad, "Yuc\\.") ~ "31",
      str_detect(ciudad, "Zac\\.") ~ "32",
      TRUE ~ NA_character_
    )
  )

# =============================================================
#               AGREGACION Y EXPORTACION
# =============================================================
inpc_mensual <- inpc_con_claves %>%
  filter(!is.na(cve_ent)) %>%
  group_by(cve_ent, fecha) %>%
  summarise(inpc = mean(inpc, na.rm = TRUE), .groups = "drop") %>%
  left_join(catalogo_inegi, by = "cve_ent") %>%
  select(cve_ent, entidad, fecha, inpc) %>%
  arrange(cve_ent, fecha)

glimpse(inpc_mensual)

write_csv(inpc_mensual, here("data", "processed", "inpc_ciudades.csv"))
