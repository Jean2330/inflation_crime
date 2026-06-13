# config.R

# ---- Credenciales ------------------------------------------------------------
BANXICO_TOKEN <- "fe7b674afcef76f947eec9328bcb2a5ad0369cfdd9e7361d4cddfe24e84de311"

FRED_KEY <- "5e0e3e17042c6b5c417c6ce2cda08d3e"

# ---- Periodo del panel -------------------------------------------------------
FECHA_INICIO <- "2004-01-01"
FECHA_FIN    <- "2024-12-31"

# ---- URLs de archivos crudos en GitHub ---------------------------------------
GITHUB_RAW_BASE <- "https://raw.githubusercontent.com/Jean2330/inflation_crime/main/data/raw"

URL_INPC_CRUDO  <- paste0(GITHUB_RAW_BASE, "/inpc_crudo.csv")
URL_ITAEE_CRUDO <- paste0(GITHUB_RAW_BASE, "/itaee_crudo.xlsx")
URL_CATALOGO    <- paste0(GITHUB_RAW_BASE, "/catalogo_inegi.csv")

# ---- Ruta local para outputs -------------------------------------------------
#RUTA_SALIDA <- "C:/Users/josel/OneDrive/Escritorio/CIDE/4to semestre/Macroeconometrics/Trabajo final/datos"
RUTA_SALIDA <- "/Users/jeanmarroquin/Documents/Economia/4to Semestre/Macroeconometria/crimen/output"

# ---- Catalogo de entidades INEGI (fuente unica, compartida por todos) --------
catalogo_inegi <- tibble::tibble(
  cve_ent = stringr::str_pad(1:32, width = 2, pad = "0"),
  entidad = c("Aguascalientes", "Baja California", "Baja California Sur",
              "Campeche", "Coahuila de Zaragoza", "Colima", "Chiapas",
              "Chihuahua", "Ciudad de México", "Durango", "Guanajuato",
              "Guerrero", "Hidalgo", "Jalisco", "México", "Michoacán de Ocampo",
              "Morelos", "Nayarit", "Nuevo León", "Oaxaca", "Puebla",
              "Querétaro", "Quintana Roo", "San Luis Potosí", "Sinaloa",
              "Sonora", "Tabasco", "Tamaulipas", "Tlaxcala",
              "Veracruz de Ignacio de la Llave", "Yucatán", "Zacatecas")
)

# ---- Series de expectativas de inflacion (12 horizontes mensuales) -----------
# La suma de estas 12 series produce la inflacion anual esperada acumulada.
# na.rm = FALSE en rowSums garantiza NA cuando algun horizonte no esta disponible.
SERIES_EXPECTATIVAS <- c(
  "SR14230",  # Mediana expectativa de inflacion mensual, horizonte t+1
  "SR14237",  # Mediana expectativa de inflacion mensual, horizonte t+2
  "SR14244",  # Mediana expectativa de inflacion mensual, horizonte t+3
  "SR14251",  # Mediana expectativa de inflacion mensual, horizonte t+4
  "SR14258",  # Mediana expectativa de inflacion mensual, horizonte t+5
  "SR14265",  # Mediana expectativa de inflacion mensual, horizonte t+6
  "SR14272",  # Mediana expectativa de inflacion mensual, horizonte t+7
  "SR14279",  # Mediana expectativa de inflacion mensual, horizonte t+8
  "SR14286",  # Mediana expectativa de inflacion mensual, horizonte t+9
  "SR14293",  # Mediana expectativa de inflacion mensual, horizonte t+10
  "SR14300",  # Mediana expectativa de inflacion mensual, horizonte t+11
  "SR14307"   # Mediana expectativa de inflacion mensual, horizonte t+12
)