rm(list = ls())
cat("\014")

# =============================================================
# LIBRERIAS
# =============================================================
library(readxl)
library(tidyverse)
library(lubridate)
library(zoo)
library(here)

# =============================================================
# CREDENCIALES
# =============================================================
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "config.R"))

# =============================================================
# LECTURA
# =============================================================
# read_excel no acepta URLs, se descarga a un archivo temporal primero
archivo_temporal <- tempfile(fileext = ".xlsx")
download.file(URL_ITAEE_CRUDO, destfile = archivo_temporal, mode = "wb", quiet = TRUE)
itaee_raw <- read_excel(archivo_temporal, col_names = FALSE)

# Busqueda automatica de encabezados
fila_anios <- which(apply(itaee_raw, 1, function(x) any(grepl("2004", as.character(x)))))[1]

# Detiene la ejecucion con mensaje claro si el Excel no tiene el formato esperado
stopifnot("No se encontro la fila con el anio 2004. Revisa el formato del archivo." = !is.na(fila_anios))

anios     <- na.locf(unlist(itaee_raw[fila_anios, -1]))
trimestres <- unlist(itaee_raw[fila_anios + 1, -1])
nombres_col <- make.unique(c("entidad", paste0(anios, "_T", trimestres)))

# =============================================================
#                     TRANSFORMACION
# =============================================================
itaee <- itaee_raw[-c(1:(fila_anios + 1)), ]
names(itaee) <- nombres_col

itaee_largo <- itaee %>%
  filter(!is.na(entidad), !str_detect(entidad, "Nota|Fuente|nacional|Nacional")) %>%
  pivot_longer(-entidad, names_to = "periodo", values_to = "itaee") %>%
  filter(str_detect(periodo, "^\\d{4}")) %>%
  mutate(
    itaee    = suppressWarnings(as.numeric(itaee)),
    anio     = as.integer(str_extract(periodo, "^\\d{4}")),
    trim_str = str_trim(str_remove(periodo, "^.*_T")),
    trim = case_when(
      trim_str %in% c("I",   "1", "01") ~ 1,
      trim_str %in% c("II",  "2", "02") ~ 2,
      trim_str %in% c("III", "3", "03") ~ 3,
      trim_str %in% c("IV",  "4", "04") ~ 4,
      TRUE ~ NA_real_
    ),
    fecha = make_date(year = anio, month = (trim - 1) * 3 + 1, day = 1)
  ) %>%
  drop_na(fecha) %>%
  select(entidad, fecha, itaee) %>%
  arrange(entidad, fecha)

# Interpolacion lineal trimestral a mensual
itaee_mensual <- itaee_largo %>%
  group_by(entidad) %>%
  complete(fecha = seq.Date(min(fecha), max(fecha), by = "month")) %>%
  mutate(itaee = approx(x    = fecha[!is.na(itaee)],
                        y    = itaee[!is.na(itaee)],
                        xout = fecha,
                        rule = 2)$y) %>%
  ungroup()

# =============================================================
#               CRUCE CON CATALOGO Y EXPORTACION
# =============================================================
itaee_mensual <- itaee_mensual %>%
  left_join(catalogo_inegi, by = "entidad") %>%
  select(cve_ent, entidad, fecha, itaee) %>%
  filter(!is.na(cve_ent),
         fecha >= "2004-01-01",
         fecha <= "2024-12-01")

glimpse(itaee_mensual)

write_csv(itaee_mensual, here("data", "processed", "itaee.csv"))
