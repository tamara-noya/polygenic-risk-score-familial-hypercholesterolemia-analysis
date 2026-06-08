######################################################################################################################
# TFM: Impacto del Riesgo Poligénico en Hipercolesterolemia Familiar
# Autora:     Tamara Noya Mosquera
# Directores: Dra. Teresa Padró - Institut de Recerca Sant Pau (IR Sant Pau)
#             Ariel Ernesto Cagariaga - Universidad Internacional de la Rioja (UNIR)
# Programa:   Máster en Bioinformática, UNIR (curso 2025-2026)
#
# Dependencias: readxl, dplyr, tidyr, anytime, gtsummary, gt, ggplot2, ggsignif, lme4, performance, pROC, tibble, 
#               survival, survminer, coxme, scales, patchwork
######################################################################################################################


######################################################################################################################
# BLOQUE 1 — CONFIGURACIÓN, CARGA Y PREPARACIÓN DE DATOS
######################################################################################################################

# ==============================================================================
# 1.1. ENTORNO DE TRABAJO
# ==============================================================================

setwd("C:/Users/LG/Documents/BIOINFO UNIR/TFM/Analisis_R/Datos")

# Directorio de salida para tablas y figuras 
ruta_out <- "C:/Users/LG/Documents/BIOINFO UNIR/TFM/Analisis_R/Figuras/Figuras_R_TMF"

# Semilla de reproducibilidad — garantiza resultados idénticos en cada ejecución.
# Referencia en Métodos: "Los análisis se realizaron con R (v4.x); semilla set.seed(2026)."
set.seed(2026)

# ==============================================================================
# 1.2. LIBRERÍAS
# ==============================================================================

library(readxl)      # Carga de archivos Excel 
library(dplyr)       # Manipulación de datos 
library(tidyr)       # Transformación de datos 
library(anytime)     # Conversión automática de strings a fechas 
library(gtsummary)   # Tablas descriptivas y de regresión con formato publicable
library(gt)          # Formato avanzado y exportación de tablas a HTML
library(ggplot2)     # Sistema de gráficos por capas 
library(ggsignif)    # Brackets de significancia estadística en gráficos ggplot
library(lme4)        # Modelos lineales/generalizados mixtos (glmer — GLMM)
library(performance) # Diagnóstico de modelos mixtos (VIF, colinealidad)
library(pROC)        # Curvas ROC, AUC y comparación de modelos (DeLong)
library(tibble)      # rownames_to_column y manipulación de tibbles
library(survival)    # Análisis de supervivencia: Surv(), coxph(), survfit()
library(survminer)   # Visualización de curvas Kaplan-Meier (ggsurvplot)
library(coxme)       # Modelos de Cox con frailty log-normal (estructura familiar)
library(scales)      # Formatos de ejes en ggplot 
library(patchwork)   # Composición de múltiples gráficos ggplot en una figura

# Tema compacto para tablas gtsummary (reduce el espaciado por defecto)
theme_gtsummary_compact()

# ==============================================================================
# 1.3. CONSTANTES DE ESTILO
# ==============================================================================
# Centralizar colores y temas garantiza consistencia visual en todo el script.

# ── Paletas de colores ─────────────────────────────────────────────────────────
# Grupo principal (FH vs No FH)
colores_fh    <- c("No FH" = "#B2BABB", "FH" = "#1F5299")
# Grupo de estudio (3 grupos: No FH / FH Sin Evento / FH Con Evento)
colores_grupo <- c("No FH" = "#B2BABB", "FH No Evento" = "#8FA4B5", "FH Evento" = "#3182CE")
# Subgrupos FH por evento (para gráficos dentro de la cohorte FH)
colores_ecv   <- c("FH No Evento" = "#8FA4B5", "FH Evento" = "#3182CE")
# Fenotipo de ECV (Sin Evento / Precoz <65a / Tardío ≥65a)
colores_fenotipo <- c("Sin Evento" = "#8FA4B5", "Precoz" = "#5DADE2", "Tardío" = "#2B6CB0")
# Sexo
colores_sexo  <- c("Mujer" = "#A593B1", "Hombre" = "#7FA493")
# Quintiles del PRS (gradiente de riesgo bajo → alto)
colores_quintile <- c("Q1" = "#4CB97B", "Q2" = "#FAD390", "Q3" = "#F89F73", "Q4" = "#D35400", "Q5" = "#A83232")
# Riesgo poligénico agrupado (Q1=Bajo, Q2-Q4=Intermedio, Q5=Alto)
colores_riesgo <- c("Bajo" = "#4CB97B", "Intermedio" = "#F89F73", "Alto" = "#A83232")
# Q5 vs resto (análisis binario de riesgo alto)
colores_q5  <- c("Q1-Q4" = "#80CFA3", "Q5" = "#A83232")
# VHR (Very High Risk: subconjunto de Q5 con criterios adicionales)
colores_vhr <- c("No" = "#80CFA3", "Si" = "#8B1E1E")
# SNP rs10455872 (LPA): portador del alelo G de riesgo
colores_lpa <- c("No portador" = "#EBDEF0", "Portador G" = "#7D3C98")
# Tipo de mutación causal de HF
colores_mut <- c("LDLR-Nulo" = "#2C5282", "LDLR-Defectivo" = "#4A7BB0", "LDLR-Desconocido" = "#A0AEC0", "ApoB" = "#63B3ED")
# Genotipos de SNPs (dosage 0/1/2)
colores_snp <- c("Hom. referencia (0)" = "#ECEFF1", "Heterocigoto (1)" = "#90A4AE", "Hom. riesgo (2)" = "#455A64")

# ── Variables y funciones auxiliares para figuras ───────────────────────────────────────────────────────────
# Se combina con theme() adicional en cada gráfico para ajustes específicos.
tema_base <- theme_bw(base_size = 10) +
  theme(strip.background      = element_rect(fill = "white"),
        strip.text            = element_text(face = "bold"),
        plot.subtitle         = element_text(size = 9, color = "gray40"),
        plot.caption          = element_text(size = 8,   color = "gray50"),
        panel.grid.major.x    = element_blank())

# Guarda un gráfico ggplot como PNG en ruta_out con parámetros estandarizados.
# Uso: guardar_figura(mi_plot, "nombre.png"). Los argumentos ancho/alto/resolucion tienen defaults pero se pueden sobreescribir.
guardar_figura <- function(plot, nombre_archivo, ancho = 7, alto = 4.5, resolucion = 300) {
  ggsave(filename = file.path(ruta_out, nombre_archivo),
         plot     = plot,
         width    = ancho,
         height   = alto,
         dpi      = resolucion)
  cat("✓ Guardada:", file.path(ruta_out, nombre_archivo), "\n")
}

# ── Funciones auxiliares para tablas gt ───────────────────────────────────────
# Aplica estilo tipográfico unificado a un objeto gt.
# Uso: datos %>% gt() %>% gt_estilo("**Título**", "*Subtítulo opcional*")
gt_estilo <- function(gt_obj, titulo, subtitulo = NULL) {
  gt_obj %>%
    tab_header(
      title    = md(titulo),
      subtitle = if (!is.null(subtitulo)) md(subtitulo) else NULL
    ) %>%
    tab_options(
      table.font.size           = px(11),             
      data_row.padding          = px(3),               
      heading.align             = "center",
      heading.padding           = px(2),
      column_labels.font.weight = "bold",
      footnotes.padding         = px(2),
      footnotes.font.size       = px(10)
    )
}

# Guarda un objeto gt como archivo HTML en ruta_out.
# Uso: gt_obj %>% guardar_gt("NombreTabla.html")
guardar_gt <- function(gt_obj, nombre_archivo) {
  gtsave(gt_obj, filename = nombre_archivo, path = ruta_out)
  cat("✓ Guardada:", file.path(ruta_out, nombre_archivo), "\n")
}

# ==============================================================================
# 1.4. CARGA DEL DATASET
# ==============================================================================

datos_raw <- read_excel("20260420_Datos_FH_PRS.xlsx")

# ==============================================================================
# 1.5. LIMPIEZA Y TRANSFORMACIÓN DE VARIABLES
# ==============================================================================

# Convierte strings con coma decimal a numérico 
limpiar_comas <- function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
}

# Se trabaja sobre una copia (datos) para preservar el raw original intacto.
# NOTA: el warning "NAs introduced by coercion" en across() de fechas es esperado
# y está controlado por F_muerte_desconocida. No indica pérdida de datos relevantes.
datos <- datos_raw %>%
  mutate(
    
    # ── A. Variables bioquímicas y de tratamiento ─────────────────────────────
    # Decimales con coma en la fuente → convertir a numérico estándar
    across(
      c(GRS, CT_0, TG_0, cHDL_0, cLDL_0, ApoA_0, ApoB_0, LpA_0,
        TSH_0, PCR_0, Glucosa_0,
        CT_1, TG_1, cHDL_1, cLDL_1, Glucosa_1, Creatinina_1,
        AñosTtoEstatinas_0, AñosEzetim_0,
        AñosTtoEstatinas_1, AñosEzetim_1, AñosTtoPCSK9_1),
      limpiar_comas),
    
    # ── B. Fechas de muerte desconocidas ──────────────────────────────────────
    # Algunos registros tienen "?" en F_muerte o Edad_muerte (fecha desconocida).
    # Se capturan antes de la conversión numérica para no perder esos casos.
    F_muerte_desconocida = factor(
      case_when(
        as.character(F_muerte)    == "?" |
          as.character(Edad_muerte) == "?"        ~ "Si",
        Muerte == "Si" & is.na(F_muerte) &
          is.na(Edad_muerte)                      ~ "Si",
        TRUE                                      ~ "No"),
      levels = c("No", "Si")),
    
    # ── C. Variables de edad: asegurar numérico ───────────────────────────────
    across(c(Edad_2025, Edad_ECV, Edad_muerte, Edad_baja, Edad_inclusion),
           ~ as.numeric(as.character(.x))),
    
    # ── D. Variables de fecha: parsear con anydate ────────────────────────────
    # anydate() tolera múltiples formatos de fecha y devuelve NA para "?"
    across(c(F_nacimiento, F_inclusion, F_ECV, F_muerte, F_baja,
             Fecha_0, Fecha_1, InicioTtoEstat_0, InicioEzetim_0,
             InicioTtoEstat_1, InicioEzetim_1, InicioTtoPCSK9_1),
           ~ anydate(.x)),
    
    # ── E. Variables categóricas ordenadas ───────────────────────────────────
    # Los niveles definen el orden en tablas y gráficos
    Sexo        = factor(Sexo,        levels = c("Mujer", "Hombre")),
    Grupo       = factor(Grupo,       levels = c("No FH", "FH")),
    Momento_ECV = factor(Momento_ECV, levels = c("preinclusion", "postinclusion")),
    Quintile    = factor(Quintile,    levels = c("Q1", "Q2", "Q3", "Q4", "Q5")),
    
    # VHR: normalizar "YES"/"NO" desde la fuente original
    VHR = factor(case_when(trimws(VHR) == "YES" ~ "Si", trimws(VHR) == "NO"  ~ "No", TRUE ~ NA_character_),
                 levels = c("No", "Si")),
    
    # HTA/DM: colapsar diagnóstico pre/post-inclusión en binario (Si/No)
    HTA_bin = factor( case_when(HTA %in% c("Diag preincl", "Diag postincl") ~ "Si", HTA == "No" ~ "No", TRUE ~ NA_character_),
                      levels = c("No", "Si")),
    
    DM_bin = factor( case_when(DM %in% c("Diag preincl", "Diag postincl") ~ "Si", DM == "No" ~ "No", TRUE ~ NA_character_),
                     levels = c("No", "Si")),
    
    # ── F. Variables dicotómicas Si/No (múltiples columnas) ──────────────────
    # La fuente original codifica indistintamente como 0/1, "Si"/"No", "si"/"no"
    across(
      c(EstudioGenPositivo_0, MutRLDL_0, MutApoB_0, ECV, Muerte, Baja, Resiliente, HTA_preinclusion, HTA_postinclusion,
        DM_preinclusion, DM_postinclusion, TtoHipolip_0, TtoEstatinas_0, TtoEzetimiba_0, TtoPCSK9_0, TtoResinas_0, 
        TtoECV_0, AAS_0, TtoHipolip_1, TtoEstatinas_1, TtoEzetimiba_1, TtoPCSK9_1, TtoResinas_1, TtoECV_1, AAS_1),
      ~ factor(case_when(.x %in% c(1, "Si", "si") ~ "Si",
                         .x %in% c(0, "No", "no") ~ "No",
                         TRUE ~ NA_character_),
               levels = c("No", "Si"))),
    
    # ── G. Variables categóricas multinivel ───────────────────────────────────
    ID_Familia       = factor(ID_Familia),
    Cod_parentesco   = factor(Cod_parentesco),
    Parentesco       = factor(Parentesco),
    Orden_parentesco = factor(Orden_parentesco),
    CodMutRLDL_0     = factor(CodMutRLDL_0),
    TipAleloRLDL_0   = factor(TipAleloRLDL_0),
    CodMutAPOB_0     = factor(CodMutAPOB_0),
    TipAleloAPOB_0   = factor(TipAleloAPOB_0),
    Tipo_ECV         = factor(Tipo_ECV),
    Causa_muerte     = factor(Causa_muerte),
    HTA              = factor(HTA),
    DM               = factor(DM),
    
    # ── H. Variables derivadas nuevas ─────────────────────────────────────────
    
    # H.1. Edad observada: última edad conocida del paciente
    #      Fallecidos → edad al fallecimiento; bajas → edad a la baja; activos → edad calculada a fecha de corte (2025)
    Edad_obs = case_when(
      Muerte == "Si" ~ Edad_muerte,
      Baja   == "Si" ~ Edad_baja,
      TRUE           ~ Edad_2025),
    
    # H.2. Grupo de estudio: variable analítica principal de estratificación
    Grupo_estudio = factor(
      case_when(Grupo == "No FH"            ~ "No FH",
                Grupo == "FH" & ECV == "No" ~ "FH No Evento",
                Grupo == "FH" & ECV == "Si" ~ "FH Evento",
                TRUE                        ~ NA_character_),
      levels = c("No FH", "FH No Evento", "FH Evento")),
    
    # H.3. Q5 vs resto: variable binaria para análisis de alto riesgo poligénico
    Quintile_5 = factor(ifelse(Quintile == "Q5", "Si", "No"), levels = c("No", "Si")),
    
    # H.4. Riesgo poligénico agrupado (categorización clínica del PRS):
    #      Bajo = Q1 (percentil 0-20), Intermedio = Q2-Q4, Alto = Q5 (percentil 80-100)
    Riesgo_poligenico = factor(
      case_when(Quintile == "Q1"                     ~ "Bajo",
                Quintile %in% c("Q2", "Q3", "Q4")   ~ "Intermedio",
                Quintile == "Q5"                     ~ "Alto",
                TRUE                                 ~ NA_character_),
      levels = c("Bajo", "Intermedio", "Alto")),
    
    # H.5. Resiliente_calc: FH genotípico sin ECV a edad ≥65 años
    #      Edad de corte 65a: umbral clínico estándar para ECV precoz/tardío en HF
    Resiliente_calc = factor(
      case_when(
        EstudioGenPositivo_0 == "Si" & ECV == "No" & Edad_obs >= 65 ~ "Si",
        EstudioGenPositivo_0 == "Si" & (ECV == "Si" | Edad_obs < 65) ~ "No",
        TRUE ~ NA_character_),
      levels = c("No", "Si")),
    
    # H.6. Fenotipo ECV en FH genotípicos: Sin Evento / Precoz (<65a) / Tardío (≥65a)
    Fenotipo_ECV = factor(
      case_when(
        EstudioGenPositivo_0 == "Si" & ECV == "Si" & Edad_ECV < 65  ~ "Precoz",
        EstudioGenPositivo_0 == "Si" & ECV == "Si" & Edad_ECV >= 65 ~ "Tardío",
        EstudioGenPositivo_0 == "Si" & ECV == "No"                  ~ "Sin Evento",
        TRUE ~ NA_character_),
      levels = c("Sin Evento", "Precoz", "Tardío")),
    
    # H.7. ID_cluster: variable de agrupación para GLMM (efectos aleatorios)
    #      Problema: los 355 singletons (F000) comparten el mismo código de familia. Si se usa ID_Familia directamente, 
    #      lme4 los agrupa en un "superclúster" de 355 individuos → efecto aleatorio espurio y sobreestimado.
    #      Solución: los F000 reciben su propio ID individual (ID_FHF); las familias reales (F001-F304) mantienen su ID_Familia compartido.
    #      NOTA: En modelos Cox con coxme, ID_Familia funciona directamente porque coxme trata singletons correctamente de forma interna.
    ID_cluster = factor(
      ifelse(as.character(ID_Familia) == "F000",
             as.character(ID_FHF),
             as.character(ID_Familia))),
    
    # H.8. TtoEstatinas_total_0: recupera uso de estatinas de dos fuentes:
    #      TtoEstatinas_0 puede estar ausente pero InicioTtoEstat_0 no (y viceversa)
    TtoEstatinas_total_0 = factor(
      case_when(
        !is.na(TtoEstatinas_0)                            ~ as.character(TtoEstatinas_0),
        is.na(TtoEstatinas_0) & !is.na(InicioTtoEstat_0) ~ "Si",
        is.na(TtoEstatinas_0) &  is.na(InicioTtoEstat_0) ~ "No",
        TRUE ~ NA_character_),
      levels = c("No", "Si")),
    
    # H.9. LLT_intensidad_0: intensidad del tratamiento hipolipemiante en inclusión
    #      Escala ordinal de 4 niveles (Lipid-Lowering Therapy intensity):
    #      0 = Sin tratamiento
    #      1 = Monoterapia con estatinas (1 fármaco)
    #      2 = Monoterapia sin estatinas (ezetimiba o resinas solos)
    #      3 = Terapia combinada (estatinas + ≥1 fármaco adicional)
    LLT_intensidad_0 = factor(
      case_when(
        TtoHipolip_0 == "No" ~ "Sin tratamiento",
        TtoEstatinas_total_0 == "Si" &
          (TtoEzetimiba_0 == "Si" | TtoResinas_0 == "Si" | TtoPCSK9_0 == "Si") ~ "Terapia combinada",
        TtoEstatinas_total_0 == "Si"                                           ~ "Monoterapia con estatinas",
        TtoHipolip_0 == "Si" & TtoEstatinas_total_0 == "No"                    ~ "Monoterapia sin estatinas",
        TRUE ~ NA_character_),
      levels = c("Sin tratamiento", "Monoterapia con estatinas",
                 "Monoterapia sin estatinas", "Terapia combinada"))
  )

# ==============================================================================
# 1.6. SUBCONJUNTOS ANALÍTICOS
# ==============================================================================
# Crear aquí todos los subsets para centralizar las definiciones.
# droplevels() en datos_fh elimina el nivel "No FH" residual de Grupo_estudio
datos_h <- filter(datos, Sexo == "Hombre")
datos_m <- filter(datos, Sexo == "Mujer")

datos_fh  <- datos %>%
  filter(Grupo == "FH") %>%
  mutate(Grupo_estudio = droplevels(Grupo_estudio),   # N = 1.086
         # Variable dependiente binaria para modelos GLMM y Cox
         ECV_bin = as.integer(ECV == "Si"))           # 0 = Sin Evento, 1 = Con Evento

datos_fh_m <- filter(datos_fh, Sexo == "Mujer")       # N = 620
datos_fh_h <- filter(datos_fh, Sexo == "Hombre")      # N = 466

cat("Cohorte completa:  N =", nrow(datos), "\n")
cat("Cohorte FH:        N =", nrow(datos_fh), "\n")
cat("Mujeres:           N =", nrow(datos_m),
    "| Total FH:",   nrow(datos_fh_m),
    "| Sin Evento:", sum(datos_fh_m$Grupo_estudio == "FH No Evento"),
    "| Con Evento:", sum(datos_fh_m$Grupo_estudio == "FH Evento"), "\n")
cat("Hombres:           N =", nrow(datos_h),
    "| Total FH:",   nrow(datos_fh_h),
    "| Sin Evento:", sum(datos_fh_h$Grupo_estudio == "FH No Evento"),
    "| Con Evento:", sum(datos_fh_h$Grupo_estudio == "FH Evento"), "\n")

# ==============================================================================
# 1.7. VALIDACIÓN POST-LIMPIEZA
# ==============================================================================

cat("\n=== Estructura general ===\n")
cat("N total:", nrow(datos), "| Variables:", ncol(datos), "\n")
cat("FH =", sum(datos$Grupo == "FH"), "| No FH =", sum(datos$Grupo == "No FH"), "\n")
cat("Familias reales (F001–F304):", nlevels(datos$ID_Familia) - 1,
    "familias,", sum(as.character(datos$ID_Familia) != "F000"), "individuos\n")
cat("Singletons (F000):", sum(as.character(datos$ID_Familia) == "F000"),
    "| Clústeres GLMM (ID_cluster):", nlevels(datos$ID_cluster), "\n")

cat("\n=== NAs en variables clave ===\n")
vars_clave <- c("GRS", "Quintile", "VHR", "ECV", "Sexo", "Edad_inclusion", "Edad_obs", "cLDL_0", "ApoB_0", "LpA_0",
                "HTA_bin", "DM_bin", "TtoEstatinas_total_0", "AñosTtoEstatinas_0", "LLT_intensidad_0", 
                "Resiliente_calc", "Fenotipo_ECV")
na_clave <- colSums(is.na(datos[vars_clave]))
if (any(na_clave > 0)) print(na_clave[na_clave > 0]) else cat("Sin NAs en variables clave.\n")

cat("\n=== Variables genéticas ===\n")
cat("Quintiles:\n");         print(table(datos$Quintile,          useNA = "always"))
cat("Riesgo poligénico:\n"); print(table(datos$Riesgo_poligenico, useNA = "always"))
cat("Q5 × VHR (VHR debe ser subconjunto estricto de Q5):\n")
print(table(Q5 = datos$Quintile_5, VHR = datos$VHR, useNA = "always"))

cat("\n=== Variables de fenotipo ===\n")
cat("Grupo de estudio:\n");  print(table(datos$Grupo_estudio,    useNA = "always"))
cat("Fenotipo ECV (FH):\n"); print(table(datos_fh$Fenotipo_ECV, useNA = "always"))
cat("Resiliente_calc:\n");   print(table(datos_fh$Resiliente_calc, useNA = "always"))
cat("Consistencia Resiliente original vs calc:\n")
print(table(Original = datos$Resiliente, Calculada = datos$Resiliente_calc,
            useNA = "always"))

cat("\n=== Tratamiento hipolipemiante en FH (N =", nrow(datos_fh), ") ===\n")
print(table(datos_fh$LLT_intensidad_0, useNA = "always"))

cat("\n=== Variables de edad ===\n")
print(summary(datos %>% select(Edad_inclusion, Edad_2025, Edad_ECV,
                               Edad_obs, Edad_muerte)))


######################################################################################################################
# BLOQUE 2 — ANÁLISIS DESCRIPTIVO
######################################################################################################################
# Objetivo: caracterizar la cohorte y comparar grupos (FH vs No FH; Sin/Con Evento)
# Tablas generadas: Tabla 1, Tabla S1, Tabla 2
# Figura generada: Boxplot GRS × Grupo × Sexo

# ==============================================================================
# 2.0. CONSTANTES COMPARTIDAS DEL BLOQUE
# ==============================================================================
# Centralizar aquí variables, etiquetas y estadísticos reutilizados en T1/S1/T2

# Variables analíticas principales 
vars_t1 <- c("Edad_inclusion", "Sexo", "CT_0", "cLDL_0", "cHDL_0", "TG_0", "ApoB_0", "LpA_0", "HTA_bin", "DM_bin",
             "TipAleloRLDL_0", "MutApoB_0", "LLT_intensidad_0", "AñosTtoEstatinas_0",
             "GRS", "Quintile", "VHR")

# Etiquetas legibles para todas las tablas (reutilizadas en T1, S1 y T2)
etiquetas_t1 <- list(
  Edad_inclusion     ~ "Edad en la inclusión (años)",
  CT_0               ~ "Colesterol total basal (mg/dL)",
  cLDL_0             ~ "cLDL basal (mg/dL)",
  cHDL_0             ~ "cHDL basal (mg/dL)",
  TG_0               ~ "Triglicéridos basales (mg/dL)",
  ApoB_0             ~ "ApoB basal (mg/dL)",
  LpA_0              ~ "Lp(a) basal (mg/dL)",
  HTA_bin            ~ "Hipertensión arterial",
  DM_bin             ~ "Diabetes mellitus",
  TipAleloRLDL_0     ~ "Tipo de alelo LDLR",
  MutApoB_0          ~ "Mutación en APOB",
  LLT_intensidad_0   ~ "Intensidad del tratamiento hipolipemiante",
  AñosTtoEstatinas_0 ~ "Años de tratamiento con estatinas",
  GRS                ~ "Genetic Risk Score (GRS)",
  Quintile           ~ "Quintil de riesgo poligénico",
  VHR                ~ "Very High Risk (VHR)")

# Estadísticos comunes: mediana (IQR) para continuas (distribuciones sesgadas confirmadas: LpA_0, GRS, 
# AñosTtoEstatinas_0), n(%) para categóricas
estadisticos_base <- list(all_continuous()  ~ "{median} ({IQR})",
                          all_categorical() ~ "{n} ({p}%)")

# Convierte p-valor en etiqueta de significancia estándar: Usado en geom_signif() para mostrar *** ** * ns de forma consistente
sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "—",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns")
}

# ==============================================================================
# 2.1. TABLA 1: Características basales FH Sin Evento vs FH Con Evento
# ==============================================================================
# Tabla principal del análisis: compara los dos grupos FH según el outcome primario (evento cardiovascular).
# Test: Wilcoxon rank-sum (continuas), Fisher exact (categóricas).
# Justificación Fisher: algunas celdas tienen n<5 (MutApoB_0, VHR).

tabla1 <- datos_fh %>%
  tbl_summary(
    by        = Grupo_estudio,
    include   = all_of(vars_t1),
    statistic = estadisticos_base,
    digits    = list(all_continuous() ~ 1, GRS ~ 3),
    missing   = "no",     # NAs documentados en notas al pie
    label     = etiquetas_t1
  ) %>%
  add_overall(last = FALSE, col_label = "**Total FH**\nN = {N}") %>%
  add_p(
    test = list(all_continuous()  ~ "wilcox.test",
                all_categorical() ~ "fisher.test"),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  bold_p() %>%
  bold_labels() %>%
  modify_header(
    label   ~ "**Variable**",
    stat_0  ~ "**Total FH**\nN = {N}",
    stat_1  ~ "**FH Sin Evento**\nN = {n}",
    stat_2  ~ "**FH Con Evento**\nN = {n}",
    p.value ~ "**p-valor**"
  ) %>%
  modify_spanning_header(c(stat_1, stat_2) ~ "**Grupo de estudio**") %>%
  modify_caption(
    "**Tabla 1.** Características basales según historia de evento cardiovascular en HF")

tabla1 %>%
  as_gt() %>%
  tab_footnote(footnote  = "Datos faltantes: FH Sin Evento n=102 (16.7%), FH Con Evento n=53 (11.1%). 
               No todos los pacientes tienen tipificación del alelo LDLR.",
               locations = cells_body(columns = label, rows = label == "Tipo de alelo LDLR")) %>%
  tab_footnote(footnote  = "Datos faltantes: FH Sin Evento n=90 (14.8%), FH Con Evento n=54 (11.3%). 
               No todos los pacientes tienen genotipificación completa de APOB.",
               locations = cells_body(columns = label, rows = label == "Mutación en APOB")) %>%
  guardar_gt("Tabla1_Descriptiva_GrupoFH.html")

# ==============================================================================
# 2.2. TABLA S1: Características basales de la cohorte completa
# ==============================================================================
# Tabla suplementaria: añade el grupo control No FH para contextualizar el perfil clínico de los pacientes FH.
# Estructura: Total cohorte | No FH | FH Total | p(No FH vs FH) | FH Sin Evento | FH Con Evento [sin p, cubierto en T1]
# Decisión de diseño: el p-valor compara No FH vs FH (2 grupos). La comparación interna FH Sin/Con Evento está en Tabla 1.

# Variables: excluir mutaciones (No FH no tiene genotipificación LDLR/ApoB)
vars_supl_t1 <- vars_t1[!vars_t1 %in% c("MutApoB_0", "TipAleloRLDL_0")]

# Filtrar etiquetas acordemente (Filter sobre lista de fórmulas)
etiquetas_supl_t1 <- Filter(function(f) !as.character(f[[2]]) %in% c("MutApoB_0", "TipAleloRLDL_0"), etiquetas_t1)

# ── Parte A: comparación principal No FH vs FH (con p-valor) ──────────────────
tabla_s1_nofh_fh <- datos %>%
  tbl_summary(
    by        = Grupo,
    include   = all_of(vars_supl_t1),
    statistic = estadisticos_base,
    digits    = list(all_continuous() ~ 1, GRS ~ 3),
    missing   = "ifany",
    missing_text = "Desconocido",
    label     = etiquetas_supl_t1
  ) %>%
  add_overall(last = FALSE, col_label = "**Total cohorte**\nN = {N}") %>%
  add_p(
    test = list(all_continuous()  ~ "wilcox.test",
                all_categorical() ~ "fisher.test"),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  bold_p() %>%
  bold_labels() %>%
  modify_header(
    label   ~ "**Variable**",
    stat_0  ~ "**Total cohorte**\nN = {N}",
    stat_1  ~ "**No FH**\nN = {n}",
    stat_2  ~ "**FH Total**\nN = {n}",
    p.value ~ "**p (No FH vs FH)**")

# ── Parte B: subgrupos FH — descriptivo únicamente, sin p ─────────────────────
tabla_s1_fh_subgrupos <- datos_fh %>%
  tbl_summary(
    by        = Grupo_estudio,    # droplevels ya aplicado en 1.6
    include   = all_of(vars_supl_t1),
    statistic = estadisticos_base,
    digits    = list(all_continuous() ~ 1, GRS ~ 3),
    missing   = "no",
    label     = etiquetas_supl_t1
  ) %>%
  modify_header(
    stat_1 ~ "**FH Sin Evento**\nN = {n}",
    stat_2 ~ "**FH Con Evento**\nN = {n}")

# ── Merge ──────────────────────────────────────────────────────────────────────
tabla_s1 <- tbl_merge(
  tbls        = list(tabla_s1_nofh_fh, tabla_s1_fh_subgrupos),
  tab_spanner = c("**Comparación principal**",
                  "**Subgrupos FH (ver Tabla 1)**")
) %>%
  modify_caption(
    "**Tabla S1.** Características basales de la cohorte completa según grupo de estudio")

tabla_s1 %>% as_gt() %>% guardar_gt("TablaS1_Descriptiva_CohorteCompleta.html")

# ==============================================================================
# 2.3. TABLA 2: Características basales FH estratificadas por sexo
# ==============================================================================
# Evalúa si las diferencias basales entre FH Sin/Con Evento son consistentes entre sexos, o si existe dimorfismo 
# sexual en el perfil de riesgo.
# Metodología: tbl_merge de dos tbl_summary independientes (M/H) para obtener p-valores específicos por sexo.
# Tests: Wilcoxon (continuas sesgadas), Fisher (categóricas, n<5 en alguna celda)
# Datos: datos_fh_m y datos_fh_h definidos en Bloque 1.6

# Variables: mismas que Tabla 1 excepto Sexo (es el factor de estratificación)
vars_t2 <- vars_t1[vars_t1 != "Sexo"]   

# Función auxiliar: evita duplicar código idéntico para mujeres y hombres. Recibe el subconjunto ya filtrado por sexo.
crear_tabla2_sexo <- function(datos_sexo) {
  datos_sexo %>%
    tbl_summary(
      by        = Grupo_estudio,   
      include   = all_of(vars_t2),
      statistic = estadisticos_base,
      digits    = list(all_continuous() ~ 1, GRS ~ 3),
      missing   = "no",            # NAs documentados en Tabla 1
      label     = etiquetas_t1
    ) %>%
    add_p(
      test = list(all_continuous()  ~ "wilcox.test",
                  all_categorical() ~ "fisher.test"),
      pvalue_fun = ~ style_pvalue(.x, digits = 3)
    ) %>%
    bold_p() %>%
    bold_labels() %>%
    modify_header(
      label   ~ "**Variable**",
      stat_1  ~ "**FH Sin Evento**\nN = {n}",
      stat_2  ~ "**FH Con Evento**\nN = {n}",
      p.value ~ "**p-valor**")
}

# Construir tablas por sexo (datos_fh_m y datos_fh_h definidos en Bloque 1.6)
tabla2_mujeres <- crear_tabla2_sexo(datos_fh_m)
tabla2_hombres <- crear_tabla2_sexo(datos_fh_h)

# Merge
tabla2_sex <- tbl_merge(
  tbls        = list(tabla2_mujeres, tabla2_hombres),
  tab_spanner = c(paste0("**Mujeres** (N = ", nrow(datos_fh_m), ")"),
                  paste0("**Hombres** (N = ", nrow(datos_fh_h), ")"))
) %>%
  modify_caption(
    "**Tabla 2.** Características basales según evento cardiovascular, estratificado por sexo")

tabla2_sex %>% as_gt() %>% guardar_gt("Tabla2_Descriptiva_FHSexo.html")

# ==============================================================================
# 2.4. DISTRIBUCIÓN DEL GRS: COMPARACIÓN ENTRE GRUPOS
# ==============================================================================
# Tres figuras complementarias que muestran la asociación del GRS con el fenotipo cardiovascular en 
# la cohorte FH y en comparación con No FH.
#
# Estructura de análisis: Global → pairwise 3 grupos → por sexo (global 3 grupos + pairwise)
#
# Diseño:
#    - Gráfico A: visión binaria No FH vs FH 
#    - Gráfico B: 3 grupos con todos los pares 
#    - Gráfico C: estratificado por sexo (solo FH (n No FH por sexo demasiado pequeña, resultados completos en texto)
#
# Tests: Wilcoxon rank-sum con corrección de Bonferroni (3 comparaciones en B)
# Se muestran todas las comparaciones incluyendo "ns", más sólido que un KW global únicamente.

# ── 2.4.1. Tests globales ─────────────────────────────────────────────────────

# Gráfico A: No FH vs FH (Wilcoxon, 2 grupos)
p_grs_bin <- wilcox.test(GRS ~ Grupo, data = datos[!is.na(datos$GRS), ])$p.value

# Gráfico B: 3 grupos - KW global + pairwise Wilcoxon Bonferroni (3 comparaciones)
kw_grs_3g <- kruskal.test(GRS ~ Grupo_estudio, data = datos[!is.na(datos$GRS), ])
pw_grs_3g <- pairwise.wilcox.test(
  x               = datos$GRS[!is.na(datos$GRS)],
  g               = datos$Grupo_estudio[!is.na(datos$GRS)],
  p.adjust.method = "bonferroni")$p.value

#Extraer los 3 pares
p_nofh_sine <- pw_grs_3g["FH No Evento", "No FH"]
p_nofh_cone <- pw_grs_3g["FH Evento",    "No FH"]
p_sine_cone <- pw_grs_3g["FH Evento",    "FH No Evento"]

# ── 2.4.2. Tests por sexo (3 grupos incluyendo No FH) ─────────────────────────
# Misma estructura que el análisis global para comparabilidad. Se usan datos_h / datos_m (incluyen No FH)

kw_grs_3g_h <- kruskal.test(GRS ~ Grupo_estudio, data = datos_h[!is.na(datos_h$GRS), ])
kw_grs_3g_m <- kruskal.test(GRS ~ Grupo_estudio, data = datos_m[!is.na(datos_m$GRS), ])

pw_grs_3g_h <- pairwise.wilcox.test(
  x               = datos_h$GRS[!is.na(datos_h$GRS)],
  g               = datos_h$Grupo_estudio[!is.na(datos_h$GRS)],
  p.adjust.method = "bonferroni")$p.value

pw_grs_3g_m <- pairwise.wilcox.test(
  x               = datos_m$GRS[!is.na(datos_m$GRS)],
  g               = datos_m$Grupo_estudio[!is.na(datos_m$GRS)],
  p.adjust.method = "bonferroni")$p.value

# Para Gráfico C (solo FH, 2 grupos por sexo)
p_grs_m <- wilcox.test(GRS ~ Grupo_estudio, data = datos_fh_m)$p.value
p_grs_h <- wilcox.test(GRS ~ Grupo_estudio, data = datos_fh_h)$p.value

cat("=== GRS entre grupos ===\n")
cat("A. No FH vs FH (Wilcoxon):               p =", round(p_grs_bin,         3), "\n")
cat("B. KW global 3 grupos:                   p =", round(kw_grs_3g$p.value, 3), "\n")
cat("B. No FH vs FH Sin Evento (Bonferroni):  p =", round(p_nofh_sine,       3), "\n")
cat("B. No FH vs FH Con Evento (Bonferroni):  p =", round(p_nofh_cone,       3), "\n")
cat("B. FH Sin vs FH Con Evento (Bonferroni): p =", round(p_sine_cone,       3), "\n")
cat("\n--- Por sexo (3 grupos incl. No FH) ---\n")
cat("Hombres — KW global:              p =", round(kw_grs_3g_h$p.value,              3), "\n")
cat("Hombres — No FH vs FH Sin Ev:     p =", round(pw_grs_3g_h["FH No Evento","No FH"], 3), "\n")
cat("Hombres — No FH vs FH Con Ev:     p =", round(pw_grs_3g_h["FH Evento","No FH"],    3), "\n")
cat("Hombres — FH Sin vs FH Con Ev:    p =", round(pw_grs_3g_h["FH Evento","FH No Evento"], 3), "\n")
cat("Mujeres — KW global:              p =", round(kw_grs_3g_m$p.value,              3), "\n")
cat("Mujeres — No FH vs FH Sin Ev:     p =", round(pw_grs_3g_m["FH No Evento","No FH"], 3), "\n")
cat("Mujeres — No FH vs FH Con Ev:     p =", round(pw_grs_3g_m["FH Evento","No FH"],    3), "\n")
cat("Mujeres — FH Sin vs FH Con Ev:    p =", round(pw_grs_3g_m["FH Evento","FH No Evento"], 3), "\n")
cat("\n--- Dentro de FH (Gráfico C) ---\n")
cat("Mujeres FH Sin vs Con Evento:     p =", round(p_grs_m, 3), "\n")
cat("Hombres FH Sin vs Con Evento:     p =", round(p_grs_h, 3), "\n")


# ── 2.4.3. Gráfico A: No FH vs FH ────────────────────────────────────────────
boxp_grs_A <- datos %>%
  filter(!is.na(GRS)) %>%
  ggplot(aes(x = Grupo, y = GRS, fill = Grupo)) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.4, color = "grey30") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_signif(comparisons = list(c("No FH", "FH")),
              annotations = sig_label(p_grs_bin),
              textsize = 3, size = 0.3, vjust = 0.3) +
  scale_fill_manual(values = colores_fh) +
  labs(title   = "Distribución del GRS: No FH vs FH",
       subtitle = paste0("Wilcoxon p = ", round(p_grs_bin, 3)),
       caption  = paste0("No FH n=", sum(datos$Grupo == "No FH"),
                         " · FH n=", sum(datos$Grupo == "FH")),
       x = NULL, y = "Genetic Risk Score (GRS)") +
  tema_base + theme(legend.position = "none")

guardar_figura(boxp_grs_A, "Boxplot_GRS_FH_vs_NoFH.png", ancho = 3.5, alto = 4.5)

# ── 2.4.4. Gráfico B: 3 grupos con todos los pares ───────────────────────────
y_max_b <- max(datos$GRS, na.rm = TRUE) * 1.35

boxp_grs_B <- datos %>%
  filter(!is.na(GRS)) %>%
  ggplot(aes(x = Grupo_estudio, y = GRS, fill = Grupo_estudio)) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.4, color = "grey30") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_signif(
    comparisons  = list(c("No FH","FH No Evento"),
                        c("No FH","FH Evento"),
                        c("FH No Evento","FH Evento")),
    annotations  = c(sig_label(p_nofh_sine),
                     sig_label(p_nofh_cone),
                     sig_label(p_sine_cone)),
    step_increase = 0.12, textsize = 3, size = 0.3, vjust = 0.3) +
  scale_fill_manual(values = colores_grupo) +
  scale_x_discrete(labels = c("No FH"        = "No FH",
                              "FH No Evento" = "FH Sin\nEvento",
                              "FH Evento"    = "FH Con\nEvento")) +
  scale_y_continuous(limits = c(NA, y_max_b)) +
  labs(title    = "Distribución del GRS por grupo de estudio",
       subtitle = paste0("KW p = ", round(kw_grs_3g$p.value, 3),
                         "  ·  Wilcoxon Bonferroni (3 comparaciones)\n",
                         "*** p<0.001  ** p<0.01  * p<0.05  ns = no significativo"),
       caption  = paste0("No FH n=", sum(!is.na(datos$GRS[datos$Grupo == "No FH"])),
                         " · FH Sin Evento n=",
                         sum(!is.na(datos$GRS[datos$Grupo_estudio == "FH No Evento"])),
                         " · FH Con Evento n=",
                         sum(!is.na(datos$GRS[datos$Grupo_estudio == "FH Evento"]))),
       x = NULL, y = "Genetic Risk Score (GRS)") +
  tema_base + theme(legend.position = "none")

guardar_figura(boxp_grs_B, "Boxplot_GRS_3Grupos.png", ancho = 5, alto = 6)

# ── 2.4.5. Gráfico C: FH Sin/Con Evento × Sexo ────────────────────────────────
# Se muestran solo FH (datos_fh) porque la N de No FH por sexo (~45-68) es insuficiente para brackets fiables. 
# Tests completos de 3 grupos disponibles en 2.4.2.
y_max_c <- max(datos_fh$GRS, na.rm = TRUE) * 1.20

signif_c <- data.frame(
  Sexo       = factor(c("Mujer","Hombre"), levels = c("Mujer","Hombre")),
  xmin       = c(1, 1), xmax = c(2, 2),
  y_position = c(y_max_c * 0.93, y_max_c * 0.93),
  label      = c(sig_label(p_grs_m), sig_label(p_grs_h)))

boxp_grs_C <- datos_fh %>%
  filter(!is.na(GRS)) %>%
  ggplot(aes(x = Grupo_estudio, y = GRS, fill = Grupo_estudio)) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.4, color = "grey30") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_signif(data        = signif_c,
              aes(annotations = label),
              xmin        = signif_c$xmin,
              xmax        = signif_c$xmax,
              y_position  = signif_c$y_position,
              manual      = TRUE, inherit.aes = FALSE,
              textsize    = 3, vjust = 0.3, size = 0.3) +
  facet_wrap(~Sexo,
             labeller = labeller(Sexo = c(
               "Mujer"  = paste0("Mujeres (n=", nrow(datos_fh_m), ")"),
               "Hombre" = paste0("Hombres (n=", nrow(datos_fh_h), ")")))) +
  scale_fill_manual(values = colores_ecv) +
  scale_x_discrete(labels = c("FH No Evento" = "FH Sin\nEvento",
                              "FH Evento"    = "FH Con\nEvento")) +
  scale_y_continuous(limits = c(NA, y_max_c)) +
  labs(title    = "Distribución del GRS por grupo de evento, estratificado por sexo",
       subtitle = paste0("Mujeres (FH): Wilcoxon p = ", round(p_grs_m, 3),
                         "  |  Hombres (FH): Wilcoxon p = ", round(p_grs_h, 3)),
       caption  = paste0("Análisis dentro de cohorte FH. Tests de 3 grupos por sexo: ",
                         "Hombres KW p=", round(kw_grs_3g_h$p.value, 3),
                         ", Mujeres KW p=", round(kw_grs_3g_m$p.value, 3)),
       x = NULL, y = "Genetic Risk Score (GRS)") +
  tema_base +
  theme(legend.position = "none", strip.text = element_text(size = 11))

guardar_figura(boxp_grs_C, "Boxplot_GRS_GrupoEvento_Sexo.png", ancho = 5.5, alto = 4.5)

# ==============================================================================
# 2.5. GRS × FENOTIPO DE ECV (Sin Evento / ECV Precoz / ECV Tardío)
# ==============================================================================
# Evalúa si el PRS discrimina no solo la presencia de ECV sino también su precocidad.
# Fenotipo: Precoz = ECV <65a, Tardío = ECV ≥65a (definido en H.6, Bloque 1). Solo FH (fenotipo ECV solo aplica a FH genotípicos)
# Estructura: global → pairwise → por sexo (misma estructura que 2.4).

# ── 2.5.1. Tests globales ─────────────────────────────────────────────────────
kw_grs_feno <- kruskal.test(GRS ~ Fenotipo_ECV, data = datos_fh)   # KW global

pw_grs_feno <- pairwise.wilcox.test(          # Wilcoxon pareados
  x               = datos_fh$GRS[!is.na(datos_fh$Fenotipo_ECV)],
  g               = datos_fh$Fenotipo_ECV[!is.na(datos_fh$Fenotipo_ECV)],
  p.adjust.method = "bonferroni")$p.value

p_sine_precoz <- pw_grs_feno["Precoz", "Sin Evento"]
p_sine_tardio <- pw_grs_feno["Tardío", "Sin Evento"]
p_prec_tardio <- pw_grs_feno["Tardío", "Precoz"]

# ── 2.5.2. Tests por sexo ─────────────────────────────────────────────────────
kw_grs_feno_m <- kruskal.test(GRS ~ Fenotipo_ECV, data = datos_fh_m)
kw_grs_feno_h <- kruskal.test(GRS ~ Fenotipo_ECV, data = datos_fh_h)

pw_grs_feno_m <- pairwise.wilcox.test(
  x = datos_fh_m$GRS[!is.na(datos_fh_m$Fenotipo_ECV)],
  g = datos_fh_m$Fenotipo_ECV[!is.na(datos_fh_m$Fenotipo_ECV)],
  p.adjust.method = "bonferroni")$p.value

pw_grs_feno_h <- pairwise.wilcox.test(
  x = datos_fh_h$GRS[!is.na(datos_fh_h$Fenotipo_ECV)],
  g = datos_fh_h$Fenotipo_ECV[!is.na(datos_fh_h$Fenotipo_ECV)],
  p.adjust.method = "bonferroni")$p.value

cat("=== GRS × Fenotipo ECV ===\n")
cat("Global - KW:               p =", round(kw_grs_feno$p.value,   3), "\n")
cat("Global - Sin Ev vs Precoz: p =", round(p_sine_precoz,         3), "(Bonferroni)\n")
cat("Global - Sin Ev vs Tardío: p =", round(p_sine_tardio,         3), "(Bonferroni)\n")
cat("Global - Precoz vs Tardío: p =", round(p_prec_tardio,         3), "(Bonferroni)\n")
cat("Mujeres - KW:              p =", round(kw_grs_feno_m$p.value, 3), "\n")
cat("Mujeres - Sin vs Precoz:   p =", round(pw_grs_feno_m["Precoz","Sin Evento"], 3), "(Bonferroni)\n")
cat("Mujeres - Sin vs Tardío:   p =", round(pw_grs_feno_m["Tardío","Sin Evento"], 3), "(Bonferroni)\n")
cat("Mujeres - Precoz vs Tardío:p =", round(pw_grs_feno_m["Tardío","Precoz"],    3), "(Bonferroni)\n")
cat("Hombres - KW:              p =", round(kw_grs_feno_h$p.value, 3), "\n")
cat("Hombres - Sin vs Precoz:   p =", round(pw_grs_feno_h["Precoz","Sin Evento"], 3), "(Bonferroni)\n")
cat("Hombres - Sin vs Tardío:   p =", round(pw_grs_feno_h["Tardío","Sin Evento"], 3), "(Bonferroni)\n")
cat("Hombres - Precoz vs Tardío:p =", round(pw_grs_feno_h["Tardío","Precoz"],    3), "(Bonferroni)\n")

# ── 2.5.3. Gráfico global ─────────────────────────────────────────────────────
y_max_feno <- max(datos_fh$GRS, na.rm = TRUE) * 1.35
n_sine   <- sum(datos_fh$Fenotipo_ECV == "Sin Evento", na.rm = TRUE)
n_precoz <- sum(datos_fh$Fenotipo_ECV == "Precoz",     na.rm = TRUE)
n_tardio <- sum(datos_fh$Fenotipo_ECV == "Tardío",     na.rm = TRUE)

boxp_grs_fenotip <- datos_fh %>%
  filter(!is.na(GRS), !is.na(Fenotipo_ECV)) %>%
  ggplot(aes(x = Fenotipo_ECV, y = GRS, fill = Fenotipo_ECV)) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.4, color = "grey30") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_signif(
    comparisons  = list(c("Sin Evento","Precoz"),
                        c("Sin Evento","Tardío"),
                        c("Precoz","Tardío")),
    annotations  = c(sig_label(p_sine_precoz),
                     sig_label(p_sine_tardio),
                     sig_label(p_prec_tardio)),
    step_increase = 0.12, textsize = 3, vjust = 0.3, size = 0.3) +
  scale_fill_manual(values = colores_fenotipo) +
  scale_x_discrete(labels = c("Sin Evento" = "Sin\nEvento",
                              "Precoz"     = "ECV Precoz\n(<65a)",
                              "Tardío"     = "ECV Tardío\n(≥65a)")) +
  scale_y_continuous(limits = c(NA, y_max_feno)) +
  labs(title    = "Distribución del GRS según fenotipo de ECV",
       subtitle = paste0("KW p = ", round(kw_grs_feno$p.value, 3),
                         "  ·  Wilcoxon Bonferroni · ns = no significativo"),
       caption  = paste0("Sin Evento n=", n_sine,
                         " · ECV Precoz n=", n_precoz,
                         " · ECV Tardío n=", n_tardio),
       x = "Fenotipo ECV", y = "Genetic Risk Score (GRS)") +
  tema_base + theme(legend.position = "none")

guardar_figura(boxp_grs_fenotip, "Boxplot_GRS_FenotipoECV.png", ancho = 4.5, alto = 5)

# ==============================================================================
# 2.6. DISTRIBUCIÓN DEL RIESGO POLIGÉNICO
# ==============================================================================
# Complementa el análisis continuo del GRS (2.4) con la distribución de las representaciones categóricas del PRS 
# entre grupos de estudio y entre sexos.
#
# Variables analizadas (ambas representan el mismo score):
#   - Riesgo_poligenico (Bajo/Intermedio/Alto): variable primaria para gráficos (3 categorías más legibles y con 
#     mayor potencia estadística al agrupar Q2-Q4)
#   - Quintile (Q1-Q5): variable secundaria para coherencia con Tabla 1 y Tabla 2 
#    
# Estructura de análisis (idéntica a 2.4 para comparabilidad):
#   Global → pairwise 3 pares → por sexo (global + pairwise)
#
# Comparaciones pairwise (3 pares con justificación científica propia):
#   · No FH vs FH Sin Evento: ¿PRS enriquecido en portadores de mutación sin ECV?
#   · No FH vs FH Con Evento: ¿PRS enriquecido en portadores con ECV?
#   · FH Sin vs FH Con Evento: ¿PRS predice ECV dentro de FH? (pregunta principal)
#
# Función chi_analisis(): ejecuta los 4 tests (global + 3 pairwise) para cualquier variable y dataset, 
# evitando repetición de código.
# Tests: Chi-cuadrado con Monte Carlo B=2000 (evita supuesto de frecuencias ≥5).

# ── 2.6.1. Funciones auxiliares para tests estadísticos ─────────────────────────────────────────────────────
# Función auxiliar: ejecuta chi² global + 3 pairwise para variable × Grupo
# Uso: chi_analisis(datos, "Riesgo_poligenico") → lista con $global, $nofh_sine, $nofh_cone, $sine_cone. 
#      El sine_cone filtra automáticamente a solo FH
chi_analisis <- function(df, var) {
  df_fh <- df %>%
    filter(Grupo_estudio != "No FH") %>%
    mutate(Grupo_estudio = droplevels(Grupo_estudio))
  list(
    global    = chisq.test(table(df[[var]], df$Grupo_estudio), simulate.p.value = TRUE, B = 2000),
    nofh_sine = chisq.test(table(droplevels(df[[var]][df$Grupo_estudio %in% c("No FH","FH No Evento")]),
                                 droplevels(df$Grupo_estudio[df$Grupo_estudio %in% c("No FH","FH No Evento")])), 
                           simulate.p.value = TRUE, B = 2000),
    nofh_cone = chisq.test(table(droplevels(df[[var]][df$Grupo_estudio %in% c("No FH","FH Evento")]),
                                 droplevels(df$Grupo_estudio[df$Grupo_estudio %in% c("No FH","FH Evento")])),
                           simulate.p.value = TRUE, B = 2000),
    sine_cone = chisq.test(table(df_fh[[var]], df_fh$Grupo_estudio), simulate.p.value = TRUE, B = 2000))
}

# ── Función de impresión compacta ──────────────────────────────────────────────
print_chi <- function(nombre, res) {
  cat(sprintf("%-28s  global=%5.3f  NoFH-Sin=%5.3f  NoFH-Con=%5.3f  Sin-Con=%5.3f\n",
              nombre,
              res$global$p.value,
              res$nofh_sine$p.value,
              res$nofh_cone$p.value,
              res$sine_cone$p.value))
}

# ── 2.6.2. Tests estadísticos ─────────────────────────────────────────────────
# Estructura: global (3 grupos) + pairwise (3 pares) × variable × sexo
# datos/datos_m/datos_h incluyen No FH; datos_fh_*/sine_cone filtra a solo FH

chi_r   <- chi_analisis(datos,   "Riesgo_poligenico")  # global
chi_r_m <- chi_analisis(datos_m, "Riesgo_poligenico")  # mujeres (incl. No FH)
chi_r_h <- chi_analisis(datos_h, "Riesgo_poligenico")  # hombres (incl. No FH)
chi_q   <- chi_analisis(datos,   "Quintile")
chi_q_m <- chi_analisis(datos_m, "Quintile")
chi_q_h <- chi_analisis(datos_h, "Quintile")

# Distribución PRS entre sexos dentro de FH
chi_r_sexo <- chisq.test(table(datos_fh$Riesgo_poligenico, datos_fh$Sexo), simulate.p.value = TRUE, B = 2000)
chi_q_sexo <- chisq.test(table(datos_fh$Quintile, datos_fh$Sexo), simulate.p.value = TRUE, B = 2000)

cat("=== Chi² distribución PRS × Grupo (Monte Carlo B=2000) ===\n")
cat(sprintf("%-28s  %s  %s  %s  %s\n", "", "Global", "NoFH-Sin", "NoFH-Con", "Sin-Con"))
print_chi("Riesgo global:",    chi_r)
print_chi("Riesgo mujeres:",   chi_r_m)
print_chi("Riesgo hombres:",   chi_r_h)
print_chi("Quintile global:",  chi_q)
print_chi("Quintile mujeres:", chi_q_m)
print_chi("Quintile hombres:", chi_q_h)
cat("Riesgo × Sexo (FH):   p =", round(chi_r_sexo$p.value, 3), "\n")
cat("Quintile × Sexo (FH): p =", round(chi_q_sexo$p.value, 3), "\n")

# ── 2.6.3. Función auxiliar: prepara datos para barras apiladas ──────────────
# Uso: prep_bar(datos, "var_eje_x", "var_relleno")
prep_bar <- function(df, var_x, var_fill) {
  df %>% filter(!is.na(.data[[var_x]]), !is.na(.data[[var_fill]])) %>%
    count(.data[[var_x]], .data[[var_fill]]) %>%
    group_by(.data[[var_x]]) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    rename(x_var = 1, fill_var = 2)
}

# Etiquetas de eje x comunes a todos los gráficos por grupo
labels_grupo <- c("No FH"        = "No FH",
                  "FH No Evento" = "FH Sin\nEvento",
                  "FH Evento"    = "FH Con\nEvento")

# ── 2.6.4. Gráfico A: Riesgo_poligenico × Grupo (global) ─────────────────────
bar_riesgo_grupo <- prep_bar(datos, "Grupo_estudio", "Riesgo_poligenico") %>%
  ggplot(aes(x = x_var, y = pct, fill = fill_var)) +
  geom_col(position = "stack", width = 0.5, alpha = 0.9) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.6, color = "white", fontface = "plain") +
  scale_fill_manual(values = colores_riesgo, name = "Riesgo\npoligénico") +
  scale_x_discrete(labels = labels_grupo) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Distribución del riesgo poligénico por grupo de estudio",
    subtitle = paste0("Chi² global p = ", round(chi_r$global$p.value, 3), "  ·  Todas las comparaciones pareadas NS"),
    caption  = paste0(
      "Bajo=Q1 · Intermedio=Q2-Q4 · Alto=Q5  ·  Monte Carlo B=2000\n",
      "Chi² pairwise: No FH vs FH Sin p=", round(chi_r$nofh_sine$p.value, 3),
      " · No FH vs FH Con p=",             round(chi_r$nofh_cone$p.value, 3),
      " · FH Sin vs FH Con p=",            round(chi_r$sine_cone$p.value, 3)),
    x = NULL, y = "Porcentaje (%)") +
  tema_base + theme(legend.position = "right")

guardar_figura(bar_riesgo_grupo, "BarChart_RiesgoPoligenico_GrupoEstudio.png", ancho = 5.5, alto = 4.5)

# ── 2.6.5. Gráfico B: Quintile × Grupo (versión detallada) ───────────────────
bar_quintile_grupo <- prep_bar(datos, "Grupo_estudio", "Quintile") %>%
  ggplot(aes(x = x_var, y = pct, fill = fill_var)) +
  geom_col(position = "stack", width = 0.6, alpha = 0.9) +
  geom_text(aes(label = ifelse(pct >= 8, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.6, color = "white", fontface = "plain") +
  scale_fill_manual(values = colores_quintile, name = "Quintil") +
  scale_x_discrete(labels = labels_grupo) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Distribución de quintiles de riesgo poligénico por grupo de estudio",
    subtitle = paste0("Chi² global p = ", round(chi_q$global$p.value, 3), "  ·  Monte Carlo B=2000"),
    caption  = paste0(
      "Chi² pairwise: No FH vs FH Sin p=", round(chi_q$nofh_sine$p.value, 3),
      " · No FH vs FH Con p=",             round(chi_q$nofh_cone$p.value, 3),
      " · FH Sin vs FH Con p=",            round(chi_q$sine_cone$p.value, 3)),
    x = NULL, y = "Porcentaje (%)") +
  tema_base + theme(legend.position = "right")

guardar_figura(bar_quintile_grupo, "BarChart_Quintile_GrupoEstudio.png", ancho = 6, alto = 4.5)

# ── 2.6.6. Gráfico C: Riesgo_poligenico × Grupo × Sexo ───────────────────────
# Figura principal del apartado. Subtítulo: p global de 3 grupos por sexo (incl. No FH) + p FH Sin vs FH Con.
bar_riesgo_grupo_sex <- datos %>%
  filter(!is.na(Riesgo_poligenico), !is.na(Grupo_estudio), !is.na(Sexo)) %>%
  count(Sexo, Grupo_estudio, Riesgo_poligenico) %>%
  group_by(Sexo, Grupo_estudio) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  ggplot(aes(x = Grupo_estudio, y = pct, fill = Riesgo_poligenico)) +
  geom_col(position = "stack", width = 0.5, alpha = 0.95) +
  geom_text(aes(label = ifelse(pct >= 8, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.6, color = "white", fontface = "plain") +
  facet_wrap(~Sexo, labeller = labeller(Sexo = c("Mujer"  = "Mujeres", "Hombre" = "Hombres"))) +
  scale_fill_manual(values = colores_riesgo, name = "Riesgo poligénico") +
  scale_x_discrete(labels = labels_grupo) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Distribución del riesgo poligénico por grupo y sexo",
    subtitle = paste0(
      "Mujeres: global p=", round(chi_r_m$global$p.value, 3), " · FH Sin vs Con p=", round(chi_r_m$sine_cone$p.value, 3), "\n",
      "Hombres: global p=", round(chi_r_h$global$p.value, 3), " · FH Sin vs Con p=", round(chi_r_h$sine_cone$p.value, 3)),
    caption  = paste0(
      "Bajo=Q1 · Intermedio=Q2-Q4 · Alto=Q5  ·  Monte Carlo B=2000\n",
      "En hombres: pairwise No FH vs FH Sin p=", round(chi_r_h$nofh_sine$p.value, 3),
      " · No FH vs FH Con p=", round(chi_r_h$nofh_cone$p.value, 3)),
    x = NULL, y = "Porcentaje (%)") +
  tema_base +
  theme(legend.position  = "bottom",
        strip.text       = element_text(face = "plain", size = 9),
        strip.background = element_rect(fill = "white", color = "white"),   
        panel.grid.major.y = element_line(linewidth = 0.3))

guardar_figura(bar_riesgo_grupo_sex, "BarChart_RiesgoPoligenico_Grupo_Sexo.png", ancho = 5.5, alto = 4.5)

# ── 2.6.7. Gráfico D: Riesgo_poligenico × Sexo (FH) ─────────────────────────
# Distribución directa entre hombres y mujeres FH 
bar_riesgo_sexo <- prep_bar(datos_fh, "Sexo", "Riesgo_poligenico") %>%
  ggplot(aes(x = x_var, y = pct, fill = fill_var)) +
  geom_col(position = "stack", width = 0.5, alpha = 0.9) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.6, color = "white", fontface = "plain") +
  scale_fill_manual(values = colores_riesgo, name = "Riesgo\npoligénico") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Distribución del riesgo poligénico por sexo\n - Cohorte FH",
    subtitle = paste0("Chi² p = ", round(chi_r_sexo$p.value, 3)),
    caption  = paste0("Bajo=Q1 · Intermedio=Q2-Q4 · Alto=Q5\n",
                      "Mujeres n=", nrow(datos_fh_m), " · Hombres n=", nrow(datos_fh_h)),
    x = NULL, y = "Porcentaje (%)") +
  tema_base + theme(legend.position = "right")

guardar_figura(bar_riesgo_sexo, "BarChart_RiesgoPoligenico_Sexo.png", ancho = 4, alto = 4.5)

# ==============================================================================
# 2.7. ANÁLISIS POST-HOC EN HOMBRES: ¿QUÉ QUINTIL IMPULSA EL EFECTO?
# ==============================================================================
# El chi² FH Sin vs FH Con en hombres es  significativo,pero no identifica qué quintil(es) son responsables. 
# Este análisis post-hoc descompone el efecto para interpretar el patrón de asociación.
#
# Dos aproximaciones complementarias:
#   A) Fisher Q1 vs Qi: ¿difiere cada quintil de Q1 en su proporción Sin/Con?
#   B) Binomial vs proporción global: ¿se aleja cada quintil de la proporción de Eventos esperada?
#
# NOTA: Los p-valores del binomial no están corregidos por múltiples comparaciones.

# ── 2.7.1. Proporción Sin/Con Evento por quintil en hombres ─────────────────
prop_quintile_h <- datos_fh_h %>%
  group_by(Quintile, Grupo_estudio) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Quintile) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

cat("=== Distribución Sin/Con Evento por quintil — Hombres FH ===\n")
print(prop_quintile_h)

# ── 2.7.2. Fisher pairwise: Qi vs Q1 (referencia) ──────────────────────────
# Pregunta A: ¿tiene cada quintil una proporción Sin/Con significativamente distinta a Q1? 
# Bonferroni corrige por las 4 comparaciones (Q2-Q5 vs Q1).
pairwise_q_h <- lapply(c("Q2","Q3","Q4","Q5"), function(q) {
  tab <- table(
    droplevels(datos_fh_h$Grupo_estudio[datos_fh_h$Quintile %in% c("Q1", q)]),
    droplevels(datos_fh_h$Quintile[datos_fh_h$Quintile      %in% c("Q1", q)]))
  p <- fisher.test(tab, simulate.p.value = TRUE, B = 2000)$p.value
  data.frame(Comparacion  = paste0("Q1 vs ", q),
             p_fisher     = round(p, 3),
             p_bonferroni = round(p.adjust(p, method = "bonferroni", n = 4), 3))
}) %>% do.call(rbind, .)

cat("\n=== Fisher pairwise Qi vs Q1 — Hombres FH ===\n")
print(pairwise_q_h)

# ── 2.7.3. Binomial vs proporción global esperada ───────────────────────────
# Pregunta B: ¿tiene cada quintil más o menos eventos de los esperados dado el % global? 
# El IC 95% binomial evalúa si la desviación es significativa.
p_eventos_h <- sum(datos_fh_h$Grupo_estudio == "FH Evento") / nrow(datos_fh_h)

binom_quintile_h <- datos_fh_h %>%
  group_by(Quintile) %>%
  summarise(
    n_total = n(),
    n_sine  = sum(Grupo_estudio == "FH No Evento"),
    n_cone  = sum(Grupo_estudio == "FH Evento"),
    pct_cone = round(n_cone / n_total * 100, 1),
    .groups = "drop") %>%
  rowwise() %>%
  mutate(
    ci_low   = round(binom.test(n_cone, n_total)$conf.int[1] * 100, 1),
    ci_high  = round(binom.test(n_cone, n_total)$conf.int[2] * 100, 1),
    p_binom  = round(binom.test(n_cone, n_total, p = p_eventos_h)$p.value, 3),
    desviacion = round(pct_cone - p_eventos_h * 100, 1)) %>%
  ungroup()

cat("\n=== Tasa de eventos por quintil vs esperada (", round(p_eventos_h*100,1), "%) ===\n")
print(binom_quintile_h %>% select(Quintile, n_total, pct_cone, ci_low, ci_high, desviacion, p_binom))

# ── 2.6.8.4. Gráfico: tasa de eventos por quintil con IC 95% ─────────────────
# Visualiza el patrón de la tasa de eventos. La línea discontinua marca la tasa global esperada.
bar_tasa_quintile_h <- binom_quintile_h %>%
  ggplot(aes(x = Quintile, y = pct_cone, fill = Quintile)) +
  geom_col(alpha = 0.85, width = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.4, color = "grey30") +
  geom_hline(yintercept = p_eventos_h * 100, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text",
           x     = 6,
           y     = p_eventos_h * 100 + 2.5,
           label = paste0("Media: ", round(p_eventos_h * 100, 1), "%"),
           size  = 2.5, color = "grey40", hjust = 1) +
  annotate("text",    # Marcar Q2 como significativo
           x = 2, y = binom_quintile_h$pct_cone[2] + 15,
           label = paste0("p=", binom_quintile_h$p_binom[2]),
           size = 2.8, color = "grey20") +
  geom_text(aes(label = paste0(pct_cone, "%"),
                y     = ci_high + 2),
            size = 2.7, color = "grey20") +
  scale_fill_manual(values = colores_quintile) +
  scale_y_continuous(limits = c(0, 90), labels = function(x) paste0(x, "%"), breaks = seq(0, 80, 20)) +
  labs(
    title    = "Tasa de evento cardiovascular por quintil - Hombres FH",
    subtitle = paste0("Línea discontinua = tasa global esperada (", round(p_eventos_h * 100, 1), "%)  · Barras = IC 95%\n", 
                      "Q2 significativamente por encima (p=", binom_quintile_h$p_binom[2],"); Q5 comparable a Q1 (Fisher p=",
                      pairwise_q_h$p_bonferroni[4], " Bonferroni)"),
    caption  = paste0("Análisis post-hoc exploratorio · Hombres FH N=", nrow(datos_fh_h),
                      " · p binomial sin corrección por múltiples comparaciones"),
    x = "Quintil de riesgo poligénico",
    y = "% con evento cardiovascular") +
  tema_base +
  theme(legend.position = "none")

guardar_figura(bar_tasa_quintile_h, "BarChart_TasaEventos_Quintile_Hombres.png", ancho = 5.5, alto = 5)


######################################################################################################################
# BLOQUE 3 - MODELOS GLMM (REGRESIÓN LOGÍSTICA MIXTA)
######################################################################################################################
#
# Objetivo: evaluar la asociación de las variables clínicas y del PRS con el evento cardiovascular en HF,
# controlando la estructura familiar mediante un efecto aleatorio por clúster familiar (ID_cluster).
#
# Justificación del diseño:
#   - El ICC > 0.10 (sección 3.1) confirma agrupamiento familiar significativo → GLMM obligatorio
#   - Variable dependiente: ECV_bin (0 = Sin Evento, 1 = Con Evento)
#   - Efecto aleatorio: (1 | ID_cluster): controla la correlación intrafamiliar 
#   - Optimizador: bobyqa (convergencia robusta en modelos mixtos con N moderada)
#   - Variables continuas estandarizadas (scale()): los OR son por 1 desviación estándar, comparables entre 
#     variables con escalas muy distintas
#
# Estructura del bloque:
#   3.1. Justificación del GLMM (ICC + DEFF)
#   3.2. Análisis univariante (OR crudas - Tabla 3)
#   3.3. Modelo base GLMM - cohorte FH completa (Tabla 4)
#   3.4. Modelos PRS - cohorte FH completa (Tabla 5 + Tabla 6)
#   3.5. Modelos PRS - estratificado por sexo (Tabla 7 + Tabla 8)
#   3.6. Test de interacción Sexo × PRS

# ==============================================================================
# 3.0. CONSTANTES COMPARTIDAS DEL BLOQUE
# ==============================================================================

# ── Variable dependiente ───────────────────────────────────────────────────────
# ECV_bin: binaria (0/1) requerida por glmer(family = binomial)
cat("ECV_bin - Cohorte FH:\n"); print(table(datos_fh$ECV_bin))
cat("ECV_bin - Mujeres FH:\n"); print(table(datos_fh_m$ECV_bin))
cat("ECV_bin - Hombres FH:\n"); print(table(datos_fh_h$ECV_bin))

# ── Optimizador GLMM ──────────────────────────────────────────────────────────
# bobyqa: derivative-free quadratic approximation. Más robusto que el optimizador por defecto (Nelder-Mead) 
# en modelos con efectos aleatorios y N moderada. maxfun = 2e5: suficiente para todos los modelos de este bloque
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ── Etiquetas de variables ────────────────────────────────────────────────────
# Centralizadas aquí para reutilizar en Tablas 3-8 sin redefinir en cada sección

# Variables del modelo base
var_labels_base <- c(
  "SexoHombre"             = "Sexo (Hombre vs. Mujer)",
  "scale(Edad_inclusion)"  = "Edad en la inclusión (por DE)",
  "scale(cLDL_0)"          = "cLDL basal (por DE)",
  "scale(cHDL_0)"          = "cHDL basal (por DE)",
  "scale(LpA_0)"           = "Lp(a) basal (por DE)",
  "HTA_binSi"              = "Hipertensión arterial (Sí vs. No)",
  "DM_binSi"               = "Diabetes mellitus (Sí vs. No)")

# Variables PRS (reutilizadas en Tablas 5, 7)
var_labels_prs <- c(
  "scale(GRS)"                     = "GRS continuo (por DE)",
  "QuintileQ2"                     = "Quintil 2 vs. Q1",
  "QuintileQ3"                     = "Quintil 3 vs. Q1",
  "QuintileQ4"                     = "Quintil 4 vs. Q1",
  "QuintileQ5"                     = "Quintil 5 vs. Q1",
  "Quintile_5Si"                   = "Quintil 5 vs. Q1-Q4",
  "Riesgo_poligenicoIntermedio"    = "Riesgo Intermedio vs. Bajo",
  "Riesgo_poligenicoAlto"          = "Riesgo Alto vs. Bajo",
  "VHRSi"                          = "VHR (Sí vs. No)")

# Filas a excluir al extraer solo las filas PRS de los modelos (varía según si el modelo incluye Sexo o no)
vars_excluir_global <- c("(Intercept)", "SexoHombre", "scale(Edad_inclusion)", "scale(cLDL_0)",
                         "scale(cHDL_0)", "scale(LpA_0)", "HTA_binSi", "DM_binSi")

vars_excluir_sexo <- vars_excluir_global[vars_excluir_global != "SexoHombre"]

# ── Funciones compartidas ─────────────────────────────────────────────────────
# Construye el modelo base + 5 modelos PRS para cualquier subconjunto de datos
# sin_sexo = TRUE: elimina Sexo del modelo (para análisis estratificados por sexo)
# Devuelve lista nombrada: $base, $grs, $quintile, $q5, $riesgo, $vhr
construir_modelos_prs <- function(df, sin_sexo = FALSE) {
  cov_base <- if (sin_sexo) {
    "scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  } else {
    "Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  }
  f_base <- as.formula(paste("ECV_bin ~", cov_base, "+ (1 | ID_cluster)"))
  list(
    base     = glmer(f_base,                                    data = df, family = binomial, control = ctrl),
    grs      = glmer(update(f_base, . ~ . + scale(GRS)),        data = df, family = binomial, control = ctrl),
    quintile = glmer(update(f_base, . ~ . + Quintile),          data = df, family = binomial, control = ctrl),
    q5       = glmer(update(f_base, . ~ . + Quintile_5),        data = df, family = binomial, control = ctrl),
    riesgo   = glmer(update(f_base, . ~ . + Riesgo_poligenico), data = df, family = binomial, control = ctrl),
    vhr      = glmer(update(f_base, . ~ . + VHR),               data = df, family = binomial, control = ctrl)
  )
}

# Extrae OR, IC 95% Wald y p-valor de un modelo glmer. Devuelve data.frame con filas nombradas por variable
tabla_ors <- function(modelo) {
  pvals <- summary(modelo)$coefficients[, "Pr(>|z|)"]
  or    <- exp(cbind(OR     = fixef(modelo),
                     confint(modelo, method = "Wald", parm = "beta_")))
  data.frame(
    OR      = round(or[, "OR"],     3),
    IC_2.5  = round(or[, "2.5 %"],  3),
    IC_97.5 = round(or[, "97.5 %"], 3),
    p_valor = ifelse(pvals < 0.001, "<0.001",
                     as.character(round(pvals, 3))))
}

# Extrae AUC y AIC de un modelo glmer
extraer_metricas <- function(modelo, nombre) {
  y_obs   <- model.response(model.frame(modelo))
  y_pred  <- fitted(modelo)
  auc_val <- as.numeric(auc(roc(y_obs, y_pred, quiet = TRUE)))
  data.frame(Modelo = nombre,
             AIC    = round(AIC(modelo), 1),
             AUC    = round(auc_val,     4))
}

# Extrae solo las filas PRS de un modelo (excluye intercepto y covariables base)
extraer_prs <- function(modelo, nombre_modelo, vars_excluir) {
  df  <- tabla_ors(modelo)
  prs <- df[!rownames(df) %in% vars_excluir, , drop = FALSE]
  if (nrow(prs) == 0) return(NULL)
  prs %>%
    tibble::rownames_to_column("Var_raw") %>%
    mutate(
      Modelo   = nombre_modelo,
      Variable = ifelse(Var_raw %in% names(var_labels_prs),
                        var_labels_prs[Var_raw], Var_raw),
      IC_95    = paste0("(", IC_2.5, " – ", IC_97.5, ")"),
      p_sig    = p_valor == "<0.001" |
        (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)) %>%
    select(Modelo, Variable, OR, IC_95, p_valor, p_sig)
}

# Construye la tabla comparativa AIC/AUC a partir de una lista de modelos
# Devuelve data.frame con Modelo, AIC, AUC y ΔAIC respecto al modelo base
comparar_modelos <- function(lista_modelos) {
  metricas <- mapply(extraer_metricas,
                     modelo  = lista_modelos,
                     nombre  = names(lista_modelos),
                     SIMPLIFY = FALSE) %>%
    do.call(rbind, .) %>%
    mutate(dAIC = round(AIC - AIC[1], 1),
           dAUC = round(AUC - AUC[1], 4))
  rownames(metricas) <- NULL
  metricas
}

# ==============================================================================
# 3.1. ANÁLISIS DEL EFECTO FAMILIAR — ICC Y DESIGN EFFECT
# ==============================================================================
# Antes de construir los modelos, se cuantifica el grado de agrupamiento familiar para comprobar si es necesario
# el uso de GLMM frente a regresión logística estándar (GLM)
#
# El ICC (Coeficiente de Correlación Intraclase) mide qué proporción de la variabilidad en ECV_bin es atribuible a 
# la familia de pertenencia. Umbral de decisión: ICC ≥ 0.10 → el agrupamiento familiar es suficientemente relevante 
# para requerir un modelo que lo controle (GLMM o GEE)
#
# El DEFF (Design Effect) cuantifica cuánto se inflan los errores estándar si se ignora el agrupamiento: 
# DEFF = 1 + (m̄ − 1) × ICC donde m̄ es el tamaño medio del clúster familiar

# ── 3.1.1. Descripción de la estructura familiar ──────────────────────────────
estructura_familiar <- datos_fh %>%
  group_by(ID_cluster) %>%
  summarise(n_miembros = n(), .groups = "drop") %>%
  summarise(
    n_clusters      = n(),
    n_singletons    = sum(n_miembros == 1),
    n_familias      = sum(n_miembros > 1),
    mediana_tamano  = median(n_miembros),
    media_tamano    = round(mean(n_miembros), 2),
    max_tamano      = max(n_miembros))

cat("=== Estructura familiar — Cohorte FH ===\n")
print(estructura_familiar)

# ── 3.1.2. Modelo nulo para estimar el ICC ────────────────────────────────────
# Modelo con solo intercepto + efecto aleatorio de familia.
# Permite estimar la varianza entre familias (σ²_u) sin contaminar con covariables.
modelo_nulo <- glmer(ECV_bin ~ 1 + (1 | ID_cluster),
                     data    = datos_fh,
                     family  = binomial(link = "logit"),
                     control = ctrl)

# Extraer varianza del efecto aleatorio e ICC en escala latente logística
# La varianza del término logístico es π²/3 ≈ 3.290 (distribución logística estándar)
var_u <- as.numeric(VarCorr(modelo_nulo)$ID_cluster)
ICC   <- var_u / (var_u + (pi^2 / 3))

# Tamaño medio del clúster e inflación de errores estándar
m_bar <- datos_fh %>%
  group_by(ID_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>% mean()

DEFF <- 1 + (m_bar - 1) * ICC

cat("\n=== ICC y Design Effect ===\n")
cat("Varianza efecto aleatorio (σ²_u):", round(var_u,       3), "\n")
cat("ICC estimado:                    ", round(ICC,         3), "\n")
cat("Tamaño medio de clúster (m̄):    ",  round(m_bar,       2), "\n")
cat("Design Effect (DEFF):            ", round(DEFF,        3), "\n")
cat("Inflación del error estándar:    ", round(sqrt(DEFF),  3), "\n")
cat("Conclusión: ICC =", round(ICC, 3),
    ifelse(ICC >= 0.10, "≥ 0.10 → efecto familiar significativo → GLMM obligatorio"))

# ── 3.1.3. Tabla ICC ──────────────────────────────────────────────────────────
tibble(
  Parametro = c(
    "Clústeres totales",
    "Individuos únicos (singletons)",
    "Familias con ≥2 miembros",
    "Tamaño medio del clúster (m̄)",
    "Tamaño máximo del clúster",
    "Varianza del efecto aleatorio (σ²_u)",
    "ICC (Coeficiente de Correlación Intraclase)",
    "Design Effect (DEFF)",
    "Inflación del error estándar (√DEFF)"),
  Valor = c(
    as.character(estructura_familiar$n_clusters),        
    as.character(estructura_familiar$n_singletons),        
    as.character(estructura_familiar$n_familias),        
    format(round(m_bar,       2), nsmall = 2),          
    as.character(estructura_familiar$max_tamano),        
    format(round(var_u,       3), nsmall = 3),           
    format(round(ICC,         3), nsmall = 3),           
    format(round(DEFF,        3), nsmall = 3),          
    format(round(sqrt(DEFF),  3), nsmall = 3)),         
  Interpretacion = c(
    "Unidades de agrupamiento en el modelo (familias + singletons con ID propio)",
    "F000: singletons con ID individual para evitar superclúster espurio",
    "Familias con correlación intrafamiliar real en el modelo",
    "Media de individuos por clúster",
    "Familia más numerosa en la cohorte",
    "Variabilidad entre familias en la escala latente logística",
    "ICC ≥ 0.10 → agrupamiento familiar significativo → GLMM obligatorio",
    "Factor de inflación efectiva del tamaño muestral por agrupamiento",
    "Factor multiplicativo sobre el error estándar si se ignora el clustering")
) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla ICC. Análisis de estructura familiar — Justificación del GLMM**") %>%
  cols_label(
    Parametro      = "Parámetro",
    Valor          = "Valor",
    Interpretacion = "Interpretación") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = 7, columns = "Valor")) %>%
  tab_footnote(
    footnote  = md("ICC = σ²_u / (σ²_u + π²/3)  ·  DEFF = 1 + (m̄ − 1) × ICC  ·  Umbral: ICC ≥ 0.10 → GLMM"),
    locations = cells_column_labels(columns = "Valor")) %>%
  tab_source_note(
    source_note = md(paste0(
      "Cohorte FH: N = ", nrow(datos_fh), ".  ", estructura_familiar$n_clusters, " clústeres (",
      estructura_familiar$n_singletons, " singletons + ", estructura_familiar$n_familias, " familias)"))) %>%
  guardar_gt("TablaICC_EstructuraFamiliar.html")


# ==============================================================================
# 3.2. ANÁLISIS UNIVARIANTE — OR NO AJUSTADOS (Tabla 3)
# ==============================================================================
# Objetivo: calcular los OR no ajustados de cada variable con ECV_bin de forma independiente, como paso previo a la
#selección de covariables del modelo base
#
# Metodología: tbl_uvregression() ejecuta un glm(family=binomial) independiente para cada variable y compila los 
# resultados en una tabla única. Los OR son por unidad de cada variable, sin estandarizar (diferente a los modelos
# ajustados donde las variables continuas se escalan con scale())
#
# Selección de variables para el modelo base (justificación clínica):
#   - Sexo: dimorfismo conocido en HF (ECV más precoz en hombres)
#   - Edad en la inclusión: factor de riesgo cardiovascular universal
#   - cLDL: biomarcador principal de HF 
#   - cHDL: factor protector independiente
#   - Lp(a): factor de riesgo independiente del LDL en HF
#   - HTA y DM: factores de riesgo clásicos con impacto independiente en HF
# Excluídas:
#   - CT y ApoB: alta colinealidad con cLDL 
#   - TG: factor de riesgo moderado, colinealidad con otros lípidos
#   - LLT y AñosTto: covariables de tratamiento, no de riesgo basal
#   - Mutación y tipo alelo: no relevantes en el modelo base

# ── 3.2.1. Tabla univariante ──────────────────────────────────────────────────
tabla3_uv <- datos_fh %>%
  select(ECV_bin, all_of(vars_t1)) %>%
  tbl_uvregression(
    method       = glm,
    y            = ECV_bin,
    method.args  = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun   = ~ style_pvalue(.x, digits = 3),
    label        = etiquetas_t1
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  modify_header(
    label    = md("**Variable**"),
    estimate = md("**OR sin ajustar**"),
    p.value  = md("**Valor p**")) %>%
  modify_column_merge(pattern = "{conf.low}, {conf.high}", rows = !is.na(conf.low)) %>%
  modify_header(conf.low = md("**IC 95%**")) %>%
  modify_caption("**Tabla 3.** Análisis univariante: factores de riesgo para ECV en la cohorte FH")

tabla3_uv %>% as_gt() %>%
  tab_footnote(
    footnote  = md("OR sin ajustar calculados mediante regresión logística binaria independiente para cada variable.
                    Variables continuas expresadas por unidad original (sin estandarizar)."),
    locations = cells_column_labels(columns = "estimate")) %>%
  tab_source_note(
    source_note = md(paste0("Cohorte FH: N = ", nrow(datos_fh), ". ECV_bin: 0 = Sin Evento (n=", sum(datos_fh$ECV_bin==0),
                            "), 1 = Con Evento (n=", sum(datos_fh$ECV_bin==1), ")"))) %>%
  guardar_gt("Tabla3_Descriptiva_UnivarianteFH.html")

# ==============================================================================
# 3.3. MODELO BASE GLMM — COHORTE FH COMPLETA (Tabla 4)
# ==============================================================================
# Modelo logístico mixto con las covariables clínicas seleccionadas en 3.2, controlando la estructura familiar
# mediante (1|ID_cluster).
#
# Variables seleccionadas:
#   Sexo + Edad_inclusion + cLDL_0 + cHDL_0 + LpA_0 + HTA_bin + DM_bin
#   Variables continuas estandarizadas: OR expresados por 1 DE para que sean comparables
#
# Diagnóstico del modelo:
#   - VIF: detecta colinealidad (umbral: VIF > 5 → problema relevante)
#   - AUC: capacidad discriminativa global del modelo base
#   - AIC: criterio de información para comparación con modelos PRS (sección 3.4)

# ── 3.3.1. Ajuste del modelo base ─────────────────────────────────────────────
modelos_global <- construir_modelos_prs(datos_fh, sin_sexo = FALSE)
modelo_base    <- modelos_global$base    # extraer para diagnóstico y tabla

summary(modelo_base)

# ── 3.3.2. Diagnóstico: colinealidad (VIF del modelo base) ────────────────────
vif_base <- check_collinearity(modelo_base)
cat("\n=== VIF — Modelo base ===\n")
print(vif_base)

max_vif     <- round(max(vif_base$VIF,       na.rm = TRUE), 2)  
max_vif_adj <- round(max(vif_base$SE_factor, na.rm = TRUE), 2)  # adj. VIF máximo (= √VIF)
cat("VIF máximo (estándar):  ", max_vif,     "\n")
cat("VIF máximo (ajustado):  ", max_vif_adj, "\n")
cat("Conclusión: ausencia de colinealidad relevante (umbral VIF < 5)\n")

# ── 3.3.3. AUC y AIC del modelo base ─────────────────────────────────────────
# Se guardan en objetos para reutilizar en sección 3.4 (tabla comparativa)
metricas_base <- extraer_metricas(modelo_base, "Base (sin PRS)")

cat("\n=== Métricas modelo base ===\n")
cat("AUC:", metricas_base$AUC, "\n")
cat("AIC:", metricas_base$AIC, "\n")

# ── 3.3.4. Tabla 4: ORs del modelo base ───────────────────────────────────────
tabla_ors(modelo_base) %>%
  tibble::rownames_to_column("Var_raw") %>%
  filter(Var_raw != "(Intercept)") %>%
  mutate(
    Variable = ifelse(Var_raw %in% names(var_labels_base), var_labels_base[Var_raw], Var_raw),
    IC_95    = paste0("(", IC_2.5, " – ", IC_97.5, ")"),
    p_sig    = p_valor == "<0.001" | 
      (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)) %>%
  select(Variable, OR, IC_95, p_valor, p_sig) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 4. Modelo GLMM base - Cohorte FH completa**",
    subtitulo = paste0("*Regresión logística mixta · N = ", nrow(datos_fh), " · AUC = ", metricas_base$AUC,
                       " · AIC = ", metricas_base$AIC, "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable",
    OR       = "OR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("OR ajustados (IC 95% Wald). Variables continuas estandarizadas (*z*-score): OR por 1 DE.<br>
                    Efecto aleatorio: (1 | ID_cluster). Optimizador: bobyqa."),
    locations = cells_column_labels(columns = "OR")) %>%
  tab_footnote(
    footnote  = md(paste0("Colinealidad verificada (VIF ajustado = √VIF): todos los VIF ajustados < ", max_vif_adj,
                          " (umbral: VIF < 5).")),
    locations = cells_title(groups = "title")) %>%
  guardar_gt("Tabla4_GLMM_Base_CohorteCompleta.html")

# ==============================================================================
# 3.4. MODELOS PRS - COHORTE FH COMPLETA (Tabla 5 + Tabla 6)
# ==============================================================================
# Los 5 modelos PRS ya fueron ajustados en 3.3.1 junto al modelo base dentro de construir_modelos_prs(). 
# Se extraen directamente de modelos_global.
#
# Tabla 5: OR de cada representación del PRS ajustada por el modelo base.
#           Responde: ¿asocia el PRS con ECV independientemente de los factores clínicos?
# Tabla 6: Comparación AIC/AUC entre modelo base y los 5 modelos PRS.
#           Responde: ¿mejora el PRS la capacidad predictiva del modelo base? Criterio: ΔAIC < −2 indica mejora relevante del ajuste.

# ── 3.4.1. Extraer ORs de los 5 modelos PRS ───────────────────────────────────
# Nombres para la tabla
nombres_prs <- c(
  grs      = "Base + GRS (continuo)",
  quintile = "Base + Quintiles",
  q5       = "Base + Q5 vs Q1-Q4",
  riesgo   = "Base + Riesgo poligénico",
  vhr      = "Base + VHR")

tabla5_datos <- lapply(names(nombres_prs), function(m) {
  extraer_prs(modelos_global[[m]],
              nombre_modelo = nombres_prs[m],
              vars_excluir  = vars_excluir_global)
}) %>%
  do.call(rbind, .) %>%
  mutate(Modelo = factor(Modelo, levels = nombres_prs))

cat("=== ORs PRS — Cohorte FH completa ===\n")
print(tabla5_datos %>% select(Modelo, Variable, OR, IC_95, p_valor))

# ── 3.4.2. Tabla 5: ORs del PRS ───────────────────────────────────────────────
tabla5_datos %>%
  select(Modelo, Variable, OR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Modelo") %>%
  gt_estilo(
    titulo    = "**Tabla 5. ORs del PRS - Cohorte FH completa**",
    subtitulo = paste0("*Efecto del riesgo poligénico ajustado por el modelo base · N = ", nrow(datos_fh), "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable PRS",
    OR       = "OR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold", color = "grey30"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote  = md("OR ajustados por Sexo, Edad, cLDL, cHDL, Lp(a), HTA y DM (IC 95% Wald).<br>
                    Variables continuas estandarizadas. Efecto aleatorio: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = "OR")) %>%
  guardar_gt("Tabla5_ORs_PRS_CFH.html")

# ── 3.4.3. Tabla 6: Comparación AIC/AUC ──────────────────────────────────────
# Renombrar la lista para que los nombres de la tabla sean legibles
modelos_global_named <- setNames(
  modelos_global, c("Base (sin PRS)", nombres_prs))

comparacion_global <- comparar_modelos(modelos_global_named)

cat("\n=== Comparación AIC/AUC — Cohorte FH completa ===\n")
print(comparacion_global)

# ── 3.4.4. Tabla 6: gt ────────────────────────────────────────────────────────
comparacion_global %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 6. Comparación de modelos PRS - Cohorte FH completa**",
    subtitulo = paste0("*AIC y AUC para el modelo base y las 5 representaciones del PRS · N = ", nrow(datos_fh), "*")) %>%
  cols_label(Modelo = "Modelo",
             AIC    = "AIC",
             AUC    = "AUC",
             dAIC   = "ΔAIC",
             dAUC   = "ΔAUC") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)) %>%
  tab_footnote(footnote  = md("ΔAIC = AIC_modelo − AIC_base."),
               locations = cells_column_labels(columns = "dAIC")) %>%
  tab_footnote(footnote  = md("ΔAUC = AUC_modelo − AUC_base.<br> Todos los modelos ajustados por (1 | ID_cluster)."),
               locations = cells_column_labels(columns = "dAUC")) %>%
  guardar_gt("Tabla6_ComparacionPRS_GrupoFH.html")

# ==============================================================================
# 3.5. MODELOS PRS - ESTRATIFICADO POR SEXO (Tabla 7 + Tabla 8)
# ==============================================================================
# Mismos 6 modelos (base + 5 PRS) para mujeres y hombres por separado
# Sexo excluido como covariable (sin_sexo = TRUE) ya que es el factor de estratificación
#
# El análisis descriptivo (Bloque 2) mostró un patrón diferencial del PRS en hombres. 
# Se evalúa si este dimorfismo persiste en el modelo ajustado
#
# Tabla 7: ORs del PRS estratificados por sexo (mujeres + hombres).
# Tabla 8: Comparación AIC/AUC por sexo.

# ── 3.5.1. Ajuste de modelos por sexo ─────────────────────────────────────────
modelos_m <- construir_modelos_prs(datos_fh_m, sin_sexo = TRUE)
modelos_h <- construir_modelos_prs(datos_fh_h, sin_sexo = TRUE)

# ── 3.5.2. ORs del PRS por sexo ───────────────────────────────────────────────
tabla7_datos <- bind_rows(
  # Mujeres
  lapply(names(nombres_prs), function(m) {
    extraer_prs(modelos_m[[m]],
                nombre_modelo = nombres_prs[m],
                vars_excluir  = vars_excluir_sexo) %>%
      mutate(Sexo = "Mujeres")
  }),
  # Hombres
  lapply(names(nombres_prs), function(m) {
    extraer_prs(modelos_h[[m]],
                nombre_modelo = nombres_prs[m],
                vars_excluir  = vars_excluir_sexo) %>%
      mutate(Sexo = "Hombres")
  })
) %>%
  mutate(Sexo  = factor(Sexo,  levels = c("Mujeres", "Hombres")),
         Modelo = factor(Modelo, levels = nombres_prs))

cat("=== ORs PRS por sexo ===\n")
print(tabla7_datos %>% select(Sexo, Modelo, Variable, OR, IC_95, p_valor))

# ── 3.5.3. Tabla 7: ORs PRS × sexo ───────────────────────────────────────────
tabla7_datos %>%
  select(Sexo, Modelo, Variable, OR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Sexo") %>%
  gt_estilo(
    titulo    = "**Tabla 7. ORs del PRS por sexo - Cohorte FH**",
    subtitulo = paste0("*Efecto del riesgo poligénico ajustado por el modelo base, estratificado por sexo*<br>",
                       "*Mujeres N = ", nrow(datos_fh_m), " · Hombres N = ", nrow(datos_fh_h), "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Modelo   = "Modelo PRS",
    Variable = "Variable",
    OR       = "OR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Modelo, Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(OR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("OR ajustados por Edad, cLDL, cHDL, Lp(a), HTA y DM (IC 95% Wald).<br>
                    Sexo excluido como covariable (análisis estratificado). Efecto aleatorio: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = "OR")) %>%
  guardar_gt("Tabla7_ORs_PRS_Sexo.html")

# ── 3.5.4. Comparación AIC/AUC por sexo ───────────────────────────────────────
modelos_m_named <- setNames(modelos_m, c("Base (sin PRS)", nombres_prs))
modelos_h_named <- setNames(modelos_h, c("Base (sin PRS)", nombres_prs))

comparacion_m <- comparar_modelos(modelos_m_named)
comparacion_h <- comparar_modelos(modelos_h_named)

cat("\n=== Comparación AIC/AUC — Mujeres ===\n"); print(comparacion_m)
cat("\n=== Comparación AIC/AUC — Hombres ===\n"); print(comparacion_h)

# ── 3.5.5. Tabla 8: AIC/AUC × sexo ───────────────────────────────────────────
left_join(
  comparacion_m %>% rename(AIC_M = AIC, AUC_M = AUC, dAIC_M = dAIC, dAUC_M = dAUC),
  comparacion_h %>% rename(AIC_H = AIC, AUC_H = AUC, dAIC_H = dAIC, dAUC_H = dAUC),
  by = "Modelo") %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 8. Comparación AIC/AUC por sexo - Modelos PRS**",
    subtitulo = paste0("*Efecto del riesgo poligénico sobre ECV, estratificado por sexo*<br>",
                       "*Mujeres N = ", nrow(datos_fh_m), " · Hombres N = ", nrow(datos_fh_h), "*")) %>%
  tab_spanner(label   = paste0("Mujeres (N = ", nrow(datos_fh_m), ")"),
              columns = c(AIC_M, AUC_M, dAIC_M, dAUC_M)) %>%
  tab_spanner(label   = paste0("Hombres (N = ", nrow(datos_fh_h), ")"),
              columns = c(AIC_H, AUC_H, dAIC_H, dAUC_H)) %>%
  cols_label(
    Modelo = "Modelo",
    AIC_M  = "AIC",  AUC_M = "AUC",
    dAIC_M = "ΔAIC", dAUC_M = "ΔAUC",
    AIC_H  = "AIC",  AUC_H = "AUC",
    dAIC_H = "ΔAIC", dAUC_H = "ΔAUC") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)) %>%
  tab_footnote(
    footnote  = md("ΔAIC = AIC_modelo - AIC_base por sexo. ΔAIC < −2 indica mejora relevante del ajuste."),
    locations = cells_column_labels(columns = "dAIC_M")) %>%
  guardar_gt("Tabla8_ComparacionAICAUC_Sexo.html")


# ==============================================================================
# 3.6. TEST DE INTERACCIÓN SEXO × PRS
# ==============================================================================
# Evalúar si el efecto del PRS sobre ECV difiere significativamente entre sexos mediante un Likelihood Ratio Test (LRT).
# H0: el efecto del PRS es igual en mujeres y hombres.
# LRT: compara modelo con efectos principales (Sexo + PRS) vs modelo con interacción (Sexo × PRS). 
# Preferible al test de Wald para modelos anidados.

# ── 3.6.1. Modelos con interacción ────────────────────────────────────────────
# update() hereda el data object exacto del modelo original
modelos_interaccion <- list(
  grs      = update(modelos_global$grs,      . ~ . + Sexo:scale(GRS)),
  quintile = update(modelos_global$quintile, . ~ . + Sexo:Quintile),
  q5       = update(modelos_global$q5,       . ~ . + Sexo:Quintile_5),
  riesgo   = update(modelos_global$riesgo,   . ~ . + Sexo:Riesgo_poligenico),
  vhr      = update(modelos_global$vhr,      . ~ . + Sexo:VHR))

# ── 3.6.2. LRT: efectos principales vs interacción ────────────────────────────
lrt_interaccion <- mapply(function(m_main, m_inter, nombre) {
  lrt  <- anova(m_main, m_inter, test = "LRT")
  data.frame(
    PRS   = nombre,
    Chi2  = round(lrt$Chisq[2], 3),
    gl    = lrt$Df[2],
    p_LRT = round(lrt$`Pr(>Chisq)`[2], 3),
    p_sig = lrt$`Pr(>Chisq)`[2] < 0.05)
},
m_main  = list(modelos_global$grs, modelos_global$quintile, modelos_global$q5,  modelos_global$riesgo, modelos_global$vhr),
m_inter = modelos_interaccion,
nombre  = c("GRS continuo", "Quintiles", "Quintil 5", "Riesgo poligénico", "VHR"),
SIMPLIFY = FALSE) %>%
  do.call(rbind, .)

cat("=== LRT Interacción Sexo × PRS ===\n")
print(lrt_interaccion %>% select(PRS, Chi2, gl, p_LRT))

# ── 3.6.3. Tabla de interacción ───────────────────────────────────────────────
lrt_interaccion %>%
  select(PRS, Chi2, gl, p_LRT, p_sig) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 9. Test de interacción Sexo × PRS**",
    subtitulo = "*Likelihood Ratio Test - modelo con efectos principales vs modelo con interacción*") %>%
  cols_label(
    PRS   = "Representación PRS",
    Chi2  = "Chi²",
    gl    = "gl",
    p_LRT = "Valor p") %>%
  cols_hide("p_sig") %>%
  cols_align(align = "left",  columns = PRS) %>%
  cols_align(align = "right", columns = c(Chi2, gl, p_LRT)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p_LRT, rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("LRT: -2·(logL_principal − logL_interacción) ~ χ² con gl = diferencia en parámetros.<br>
                    gl: grados de libertad del término de interacción (1 para variables binarias, k−1 para categóricas)."),
    locations = cells_column_labels(columns = Chi2)) %>%
  tab_source_note(
    source_note = md(paste0("Cohorte FH: N = ", nrow(datos_fh),
                            ". Modelos ajustados por Sexo, Edad, cLDL, cHDL, Lp(a), HTA, DM + (1|ID_cluster)"))) %>%
  guardar_gt("Tabla9_Interaccion_SexoPRS.html")



######################################################################################################################
# BLOQUE 4 — ANÁLISIS DE SUPERVIVENCIA (KAPLAN-MEIER + COX CON FRAILTY)
######################################################################################################################

# Objetivo: modelar el tiempo hasta el primer evento cardiovascular en la cohorte FH, complementando el análisis 
# transversal del GLMM (Bloque 3) con la dimensión temporal.
#
# Diseño del modelo:
#   - Variable de respuesta: Surv(t_evento, ECV_bin), donde t_evento = tiempo desde los 18 años
#   - Tiempo desde los 18 años (no desde la inclusión): captura la exposición acumulada al riesgo desde la edad biológica 
#     de inicio del riesgo cardiovascular en HF → evita distorsión por diferentes edades de entrada al registro
#   - Efecto aleatorio familiar: frailty term (1 | ID_cluster) via coxme, análogo al GLMM
#   - Sexo: incluido como strata(Sexo) porque viola el supuesto de riesgos proporcionales (verificado mediante 
#     residuos de Schoenfeld, sección 4.2) → HR para Sexo no estimado
#   - Variables continuas estandarizadas (scale()): HR comparables entre covariables

# ==============================================================================
# 4.0. CONSTANTES COMPARTIDAS DEL BLOQUE
# ==============================================================================

# ── Etiquetas de variables del modelo base Cox ────────────────────────────────
# Sexo no incluido: estratificado via strata(Sexo), no estimado como coeficiente
var_labels_base_cox <- c(
  "scale(Edad_inclusion)"  = "Edad en la inclusión (por DE)",
  "scale(cLDL_0)"          = "cLDL basal (por DE)",
  "scale(cHDL_0)"          = "cHDL basal (por DE)",
  "scale(LpA_0)"           = "Lp(a) basal (por DE)",
  "HTA_binSi"              = "Hipertensión arterial (Sí vs. No)",
  "DM_binSi"               = "Diabetes mellitus (Sí vs. No)")

# Filas a excluir al extraer solo filas PRS (idéntico para global y estratificado por sexo)
vars_excluir_cox <- c("scale(Edad_inclusion)", "scale(cLDL_0)", "scale(cHDL_0)", "scale(LpA_0)", "HTA_binSi", "DM_binSi")

# ── Funciones compartidas ─────────────────────────────────────────────────────
# Construye modelo Cox base + 5 modelos PRS con frailty familiar
# sin_sexo = FALSE (global):    incluye strata(Sexo) en la fórmula
# sin_sexo = TRUE  (por sexo):  Sexo excluido, análisis dentro de cada sexo
# Devuelve lista nombrada: $base, $grs, $quintile, $q5, $riesgo, $vhr
construir_modelos_cox <- function(df, sin_sexo = FALSE) {
  cov_base <- if (sin_sexo) {
    "scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  } else {
    "strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  }
  f_base <- as.formula(paste("Surv(t_evento, ECV_bin) ~", cov_base, "+ (1 | ID_cluster)"))
  list(
    base     = coxme(f_base,                                    data = df),
    grs      = coxme(update(f_base, . ~ . + scale(GRS)),        data = df),
    quintile = coxme(update(f_base, . ~ . + Quintile),          data = df),
    q5       = coxme(update(f_base, . ~ . + Quintile_5),        data = df),
    riesgo   = coxme(update(f_base, . ~ . + Riesgo_poligenico), data = df),
    vhr      = coxme(update(f_base, . ~ . + VHR),               data = df)
  )
}

# Extrae HR, IC 95% (Wald: ±1,96·SE) y p-valor de un modelo coxme
# Devuelve data.frame con filas nombradas por variable
extraer_hrs <- function(modelo) {
  cf <- summary(modelo)$coefficients
  data.frame(
    HR      = round(cf[, "exp(coef)"],                              3),
    IC_2.5  = round(exp(cf[, "coef"] - 1.96 * cf[, "se(coef)"]),   3),
    IC_97.5 = round(exp(cf[, "coef"] + 1.96 * cf[, "se(coef)"]),   3),
    p_valor = ifelse(cf[, "p"] < 0.001, "<0.001",
                     as.character(round(cf[, "p"], 3))))
}

# Extrae AIC de un modelo coxme para comparación de modelos
# coxme no implementa AIC() directamente → se calcula manualmente: AIC = -2·logLik + 2·(nº parámetros fijos + 1 varianza frailty)
extraer_aic_cox <- function(modelo, nombre) {
  ll  <- modelo$loglik["Integrated"]          # log-verosimilitud del modelo ajustado
  k   <- length(modelo$coefficients)          # nº parámetros fijos
  aic <- round(-2 * ll + 2 * (k + 1), 1)      # +1 por el parámetro de varianza del frailty
  data.frame(Modelo = nombre, AIC = aic)
}

# Construye tabla comparativa de AIC para una lista de modelos coxme.
# Devuelve data.frame con Modelo, AIC y ΔAIC respecto al modelo base.
comparar_modelos_cox <- function(lista_modelos) {
  mapply(extraer_aic_cox,
         modelo  = lista_modelos,
         nombre  = names(lista_modelos),
         SIMPLIFY = FALSE) %>%
    do.call(rbind, .) %>%
    mutate(dAIC = round(AIC - AIC[1], 1))  %>%
    {rownames(.) <- NULL; .}
}

# Extrae solo las filas PRS de un modelo coxme. Reutiliza var_labels_prs definido en Bloque 3 
extraer_prs_cox <- function(modelo, nombre_modelo, vars_excluir) {
  df  <- extraer_hrs(modelo)
  prs <- df[!rownames(df) %in% vars_excluir, , drop = FALSE]
  if (nrow(prs) == 0) return(NULL)
  prs %>%
    tibble::rownames_to_column("Var_raw") %>%
    mutate(Modelo   = nombre_modelo,
           Variable = ifelse(Var_raw %in% names(var_labels_prs), var_labels_prs[Var_raw], Var_raw),
           IC_95    = paste0("(", IC_2.5, " – ", IC_97.5, ")"),
           p_sig    = p_valor == "<0.001" |
             (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)) %>%
    select(Modelo, Variable, HR, IC_95, p_valor, p_sig)
}

# ==============================================================================
# 4.1. CREACIÓN DEL DATASET DE SUPERVIVENCIA
# ==============================================================================
# Se excluyen los individuos que fallecieron por causa desconocida sin evento documentado, ya que no es posible 
# determinar su tiempo de censura.Para el resto: tiempo al evento = Edad_ECV − 18 (si ECV) o Edad_obs − 18 (si censurado)
# Se restan 18 años porque la HF es una enfermedad congénita con riesgo cardiovascular acumulado desde la infancia; 
# anclar en los 18 años captura la exposición de riesgo desde la edad adulta temprana 
#
# t_post (desde la inclusión): reservado para el análisis de sensibilidad (4.8).
# Permite restringir el análisis a eventos incidentes (post-inclusión), que descarta el sesgo de left truncation 
# por eventos prevalentes.

datos_surv <- datos_fh %>%
  filter(!(Muerte == "Si" & is.na(Edad_muerte) & ECV != "Si")) %>%
  mutate(
    # Tiempo principal: desde los 18 años
    t_evento = case_when(
      ECV_bin == 1 ~ Edad_ECV       - 18,
      TRUE         ~ Edad_obs        - 18),
    # Tiempo desde inclusión: para análisis de sensibilidad (4.8)
    t_post   = case_when(
      ECV_bin == 1 ~ Edad_ECV       - Edad_inclusion,
      TRUE         ~ Edad_obs        - Edad_inclusion))

# Subsets por sexo (heredan t_evento y t_post)
datos_surv_m <- filter(datos_surv, Sexo == "Mujer")
datos_surv_h <- filter(datos_surv, Sexo == "Hombre")

# ── Verificación del dataset ───────────────────────────────────────────────────
cat("=== Dataset de supervivencia ===\n")
cat("N total:                 ", nrow(datos_surv),                          "\n")
cat("Eventos (ECV_bin = 1):   ", sum(datos_surv$ECV_bin),                   "\n")
cat("Censurados:              ", sum(datos_surv$ECV_bin == 0),              "\n")
cat("Tiempos negativos:       ", sum(datos_surv$t_evento < 0, na.rm=TRUE),  "\n")
cat("NAs en t_evento:         ", sum(is.na(datos_surv$t_evento)),           "\n")
cat("NAs en t_post:           ", sum(is.na(datos_surv$t_post)),             "\n")
cat("Mujeres N:               ", nrow(datos_surv_m),                        "\n")
cat("Hombres N:               ", nrow(datos_surv_h),                        "\n")

# ==============================================================================
# 4.2. SUPUESTO DE RIESGOS PROPORCIONALES - TEST DE SCHOENFELD
# ==============================================================================
# El modelo de Cox asume que el Hazard Ratio de cada covariable es constante a lo largo del tiempo de seguimiento 
# (supuesto PH). Si se viola, el HR estimado es un promedio ponderado de un efecto que cambia con el tiempo y lleva a una
# interpretación engañosa.
#
# Verificación mediante residuos de Schoenfeld (cox.zph):
#   - H0: no hay correlación entre los residuos y el tiempo → supuesto PH cumplido
#   - p < 0.05: supuesto PH violado para esa variable
#
# Nota: cox.zph() requiere un objeto coxph (no coxme) → se ajusta un modelo auxiliar sin frailty solo para el diagnóstico. 
# Los modelos finales (sección 4.4+) incluyen el frailty term (1 | ID_cluster).

# ── 4.2.1. Modelo auxiliar para test de Schoenfeld ────────────────────────────
cox_ph_diag <- coxph(
  Surv(t_evento, ECV_bin) ~ Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin,
  data = datos_surv)

# ── 4.2.2. Test de proporcionalidad de riesgos ────────────────────────────────
ph_test <- cox.zph(cox_ph_diag, transform = "km")

cat("=== Test de Schoenfeld — Supuesto de riesgos proporcionales ===\n")
print(ph_test)

# Gráfico de residuos de Schoenfeld: Línea horizontal en 0 → supuesto cumplido; tendencia temporal → violación
plot(ph_test, var = "Sexo",
     main = "Residuos de Schoenfeld - Sexo",
     xlab = "Tiempo (años desde 18)", ylab = "Beta(t)")
abline(h = 0, lty = 2, col = "grey50")

cat("\n=== Interpretación Schoenfeld ===\n")
cat("Sexo: p =", format(ph_test$table["Sexo", "p"], digits=3), "→ PH violado → strata(Sexo) en todos los modelos\n")
cat("Edad: p <2e-16 → violación esperada (sesgo supervivencia selectiva)\n")
cat("cLDL: p = 0.003 → violación coherente con cambio tratamiento en seguimiento\n")
cat("HTA:  p = 0.008 → borderline; mantenida como time-invariant\n")
cat("Decisión: strata(Sexo) + resto como covariables + limitación documentada\n")

# ==============================================================================
# 4.3. CURVAS DE KAPLAN-MEIER
# ==============================================================================
# Objetivo: visualizar la supervivencia libre de ECV según el nivel de riesgo poligénico, para evaluar descriptivamente 
# si el PRS discrimina el tiempo al evento en la cohorte FH.
#
# Figura 2A: KM × Riesgo_poligenico (cohorte FH global, 3 grupos). Pregunta: ¿se separan las curvas por nivel de riesgo poligénico?
# Figura 2B: KM × Quintile en hombres (5 grupos). Pregunta: ¿se observa algún patrón?
#
# Log-rank test: prueba global de igualdad de curvas de supervivencia.
# No asume proporcionalidad de riesgos → válido aunque cox.zph detectara violaciones en algunas variables del modelo ajustado.
#
# Nota: Los p-valores mostrados corresponden al log-rank test.

# ── Función para guardar ggsurvplot con risk table ────────────────────────────
# guardar_figura() usa ggsave() que no es compatible con ggsurvplot + risk.table.
# Esta función auxiliar usa png() + print() para preservar el risk table.
guardar_km <- function(km_obj, nombre_archivo, ancho = 8.5, alto = 7) {
  ruta_completa <- file.path(ruta_out, nombre_archivo)
  png(ruta_completa, width = ancho, height = alto, units = "in", res = 300)
  print(km_obj)
  dev.off()
  cat("✓ Guardada:", ruta_completa, "\n")
}

# ── Parámetros estéticos comunes a todas las figuras KM ───────────────────────
opciones_km <- list(                     # Sin CI marcados
  conf.int          = FALSE,           
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,              # p-valor log-rank dentro del gráfico
  pval.size         = 3.5,
  pval.method       = TRUE,
  pval.coord        = c(2, 0.25),      
  xlab              = "Tiempo desde los 18 años (años)",
  ylab              = "Probabilidad libre de ECV",
  ggtheme           = theme_bw(base_size = 11) +
    theme(panel.grid.minor  = element_blank(),
          legend.position   = "right",
          legend.key.size   = unit(0.4, "cm")),
  tables.theme      = theme_cleantable(),
  fontsize          = 3.2,
  risk.table.y.text = FALSE)

# Con CI marcados
opciones_km_ci <- modifyList(opciones_km, list(conf.int = TRUE, conf.int.alpha = 0.1))


# ── 4.3.1. Log-rank tests ──────────────────────────────────────────────────────
lr_riesgo        <- survdiff(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv)
lr_quintile      <- survdiff(Surv(t_evento, ECV_bin) ~ Quintile,   data = datos_surv)
lr_vhr           <- survdiff(Surv(t_evento, ECV_bin) ~ VHR,        data = datos_surv)
lr_q5            <- survdiff(Surv(t_evento, ECV_bin) ~ Quintile_5, data = datos_surv)
lr_h_quintile    <- survdiff(Surv(t_evento, ECV_bin) ~ Quintile, data = datos_surv_h)
lr_h_riesgo      <- survdiff(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv_h)
lr_m_quintile    <- survdiff(Surv(t_evento, ECV_bin) ~ Quintile,   data = datos_surv_m)
lr_m_riesgo      <- survdiff(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv_m)

# Función auxiliar para extraer p-valor del log-rank
p_lr <- function(lr_obj) {
  round(1 - pchisq(lr_obj$chisq, df = length(lr_obj$n) - 1), 3)
}

cat("\n=== Log-rank tests - resumen completo ===\n")
cat("Global FH:\n")
cat("  Riesgo_poligenico:    p =", p_lr(lr_riesgo),     "\n")
cat("  Quintile (Q1-Q5):     p =", p_lr(lr_quintile),   "\n")
cat("  Quintile_5 (binario): p =", p_lr(lr_q5),         "\n")
cat("  VHR:                  p =", p_lr(lr_vhr),         "\n")
cat("Mujeres FH:\n")
cat("  Quintile (Q1-Q5):     p =", p_lr(lr_m_quintile), "\n")
cat("  Riesgo_poligenico:    p =", p_lr(lr_m_riesgo),   "\n")
cat("Hombres FH:\n")
cat("  Quintile (Q1-Q5):     p =", p_lr(lr_h_quintile), "\n")
cat("  Riesgo_poligenico:    p =", p_lr(lr_h_riesgo),   "\n")

# ── 4.3.2. KM × Riesgo_poligenico - Cohorte FH global ──────────────
# Visualiza si las 3 categorías de riesgo poligénico se asocian con diferente tiempo al evento en la cohorte FH completa 
# (log-rank global).
fit_km_riesgo <- survfit(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv)
km_riesgo <- do.call(ggsurvplot, c(
  list(fit          = fit_km_riesgo,
       data         = datos_surv,
       palette      = unname(colores_riesgo),
       legend.title = "Riesgo poligénico",
       legend.labs  = c("Bajo (Q1)", "Intermedio (Q2-Q4)", "Alto (Q5)"),
       title        = "Kaplan-Meier por Categoría de Riesgo poligénico (N=1.051)"),
  opciones_km_ci))
guardar_km(km_riesgo, "KM_FH_Riesgo.png")

# ── 4.3.3. KM × Quintile - Cohorte FH global ──────────────────
fit_km_quintile <- survfit(Surv(t_evento, ECV_bin) ~ Quintile, data = datos_surv)
km_quintile <- do.call(ggsurvplot, c(
  list(fit          = fit_km_quintile,
       data         = datos_surv,
       palette      = unname(colores_quintile),
       legend.title = "Quintil PRS",
       legend.labs  = c("Q1", "Q2", "Q3", "Q4", "Q5"),
       title        = "Kaplan-Meier por Quintil de Riesgo Poligénico - Cohorte FH (N=1.051)"),
  opciones_km))
guardar_km(km_quintile, "KM_FH_Quintile.png")

# ── 4.3.4. KM × Quintile_5 binario - Cohorte FH global ────────
fit_km_q5 <- survfit(Surv(t_evento, ECV_bin) ~ Quintile_5, data = datos_surv)
km_q5 <- do.call(ggsurvplot, c(
  list(fit          = fit_km_q5,
       data         = datos_surv,
       palette      = unname(colores_q5),
       legend.title = "Quintil 5",
       legend.labs  = c("Q1-Q4", "Q5"),
       title        = "Kaplan-Meier por Quintil 5 binario - Cohorte FH (N=1.051)"),
  opciones_km_ci))
guardar_km(km_q5, "KM_FH_Q5.png")

# ── 4.3.5. KM × VHR - Cohorte FH global ─────────────────────
fit_km_vhr <- survfit(Surv(t_evento, ECV_bin) ~ VHR, data = datos_surv)
km_vhr <- do.call(ggsurvplot, c(
  list(fit          = fit_km_vhr,
       data         = datos_surv,
       palette      = unname(colores_vhr),
       legend.title = "VHR",
       legend.labs  = c("No", "Sí"),
       title        = "Kaplan-Meier por VHR (Very High Risk) - Cohorte FH (N=1.051)"),
  opciones_km_ci))
guardar_km(km_vhr, "KM_FH_VHR.png")

# ── 4.3.6. KM × Quintile - Hombres FH ─────────────────────────────
fit_km_h_quintile <- survfit(Surv(t_evento, ECV_bin) ~ Quintile, data = datos_surv_h)
km_h_quintile <- do.call(ggsurvplot, c(
  list(fit          = fit_km_h_quintile,
       data         = datos_surv_h,
       palette      = unname(colores_quintile),
       legend.title = "Quintil PRS",
       legend.labs  = c("Q1", "Q2", "Q3", "Q4", "Q5"), 
       title        = paste0("Kaplan-Meier por Quintil de Riesgo Poligénico - Hombres FH (n=", nrow(datos_surv_h), ")")),
  opciones_km))
guardar_km(km_h_quintile, "KM_H_Quintile.png")

# ── 4.3.7. KM × Riesgo_poligenico - Hombres FH ──────────────
fit_km_h_riesgo <- survfit(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv_h)
km_h_riesgo <- do.call(ggsurvplot, c(
  list(fit          = fit_km_h_riesgo,
       data         = datos_surv_h,
       palette      = unname(colores_riesgo),
       legend.title = "Riesgo poligénico",
       legend.labs  = c("Bajo (Q1)", "Intermedio (Q2-Q4)", "Alto (Q5)"),
       title        = paste0("Kaplan-Meier por Categoría de Riesgo poligénico - Hombres FH (n=", nrow(datos_surv_h), ")")),
  opciones_km_ci))
guardar_km(km_h_riesgo, "KM_H_Riesgo.png")

# ── 4.3.8. KM × Quintile - Mujeres FH ────────────────────────
fit_km_m_quintile <- survfit(Surv(t_evento, ECV_bin) ~ Quintile, data = datos_surv_m)
km_m_quintile <- do.call(ggsurvplot, c(
  list(fit          = fit_km_m_quintile,
       data         = datos_surv_m,
       palette      = unname(colores_quintile),
       legend.title = "Quintil PRS",
       legend.labs  = c("Q1", "Q2", "Q3", "Q4", "Q5"),
       title        = paste0("Kaplan-Meier por Quintil de Riesgo Poligénico - Mujeres FH (n=", nrow(datos_surv_m), ")")),
  opciones_km))
guardar_km(km_m_quintile, "KM_M_Quintile.png")

# ── 4.3.9. KM × Riesgo_poligenico - Mujeres FH ───────────────
fit_km_m_riesgo <- survfit(Surv(t_evento, ECV_bin) ~ Riesgo_poligenico, data = datos_surv_m)
km_m_riesgo <- do.call(ggsurvplot, c(
  list(fit          = fit_km_m_riesgo,
       data         = datos_surv_m,
       palette      = unname(colores_riesgo),
       legend.title = "Riesgo poligénico",
       legend.labs  = c("Bajo (Q1)", "Intermedio (Q2-Q4)", "Alto (Q5)"),
       title        = paste0("Kaplan-Meier por Categoría de Riesgo poligénico - Mujeres FH (n=", nrow(datos_surv_m), ")")),
  opciones_km_ci))
guardar_km(km_m_riesgo, "KM_M_Riesgo.png")


# ==============================================================================
# 4.4. MODELO COX BASE CON FRAILTY (Tabla Cox base)
# ==============================================================================
# Modelo de regresión de Cox con frailty familiar que estima el HR de cada covariable clínica sobre el tiempo al evento cardiovascular.
#
# Diseño:
#   - strata(Sexo): justificado por violación del supuesto PH (Schoenfeld p=5.8e-05) -> la tasa de riesgo base se 
#     estima por separado para hombres y mujeres, pero los HR de las covariables son comunes a ambos sexos
#   - (1 | ID_cluster): frailty gamma, controla la correlación intrafamiliar (análogo al efecto aleatorio del GLMM, 
#     pero en escala de tiempo)
#   - Variables continuas estandarizadas (scale()): HR por 1 DE, comparables
#
# Diagnóstico: varianza del frailty (σ² > 0 → agrupamiento familiar relevante)
# Comparación de modelos: ΔAIC de la log-verosimilitud integrada (C-index no reportado, ver justificación en 4.4.3)

# ── 4.4.1. Ajuste del modelo base Cox ─────────────────────────────────────────
modelos_cox_global <- construir_modelos_cox(datos_surv, sin_sexo = FALSE)
cox_base           <- modelos_cox_global$base

summary(cox_base)

# ── 4.4.2. Diagnóstico: varianza del frailty ──────────────────────────────────
var_frailty <- cox_base$vcoef[[1]]   # varianza del término de frailty
cat("\n=== Frailty - Modelo Cox base ===\n")
cat("Varianza del frailty (σ²):", round(var_frailty, 3), "\n")
cat("Interpretación: σ² > 0 → agrupamiento familiar relevante en supervivencia\n")

# ── 4.4.3. AIC del modelo base Cox (C-index no reportado) ──────────────────────────
# El C-index no se reporta para los modelos Cox porque t_evento (tiempo desde los 18 años) y Edad_inclusion son variables 
# estructuralmente dependientes: pacientes de mayor edad al entrar tienen simultáneamente t_evento más largo y LP más 
# bajo (coef Edad = -0.512), produciendo concordancias no interpretables independientemente del método de cálculo.
# La discriminación del modelo GLMM proporciona la medida de discriminación principal. Los modelos Cox 
# se comparan entre sí mediante ΔAIC.
cat("AIC modelo base (loglik integrada):\n")
print(cox_base$loglik)    

# ΔAIC = AIC_modelo_PRS − AIC_base: negativo indica mejora del ajuste.
aic_cox_base <- extraer_aic_cox(cox_base, "Base (sin PRS)")

cat("=== AIC modelo Cox base ===\n")
cat("AIC:", aic_cox_base$AIC, "\n")
cat("(Para comparación: ΔAIC < −2 indica mejora relevante del ajuste)\n")

# ── 4.4.4. Tabla Cox base ─────────────────────────────────────────────────────
extraer_hrs(cox_base) %>%
  tibble::rownames_to_column("Var_raw") %>%
  filter(Var_raw != "(Intercept)") %>%
  mutate(
    Variable = ifelse(Var_raw %in% names(var_labels_base_cox), var_labels_base_cox[Var_raw], Var_raw),
    IC_95    = paste0("(", IC_2.5, " - ", IC_97.5, ")"),
    p_sig    = p_valor == "<0.001" |
      (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)) %>%
  select(Variable, HR, IC_95, p_valor, p_sig) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 10. Modelo Cox base con frailty familiar - Cohorte FH**",
    subtitulo = paste0("*Cox proporcional · N = ",  cox_base$n[2], " · Estratificado por sexo",
                       " (σ² = ", round(var_frailty, 3), ")*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("HR ajustados (IC 95%: exp(β ± 1,96·SE)). Variables continuas estandarizadas (*z*-score).<br>
                    Sexo incluido como variable de estratificación: HR para Sexo no estimado. Frailty: ID_cluster"),
    locations = cells_column_labels(columns = "HR")) %>%
  tab_source_note(
    source_note = md(paste0( "Cohorte FH: N = ",  cox_base$n[2], ". Eventos = ", cox_base$n[1],
                             ". Censurados = ", cox_base$n[2] - cox_base$n[1]))) %>%
  guardar_gt("Tabla10_Cox_Base_CohorteCompleta.html")

# ==============================================================================
# 4.5. MODELOS PRS - COHORTE FH COMPLETA (Tabla 11 + comparación AIC)
# ==============================================================================
# Los 5 modelos PRS ya fueron ajustados en 4.4.1 dentro de construir_modelos_cox(). Se extraen directamente de modelos_cox_global.
#
# Tabla 11: HR de cada representación del PRS ajustada por el modelo base Cox.
#            Pregunta: ¿asocia el PRS con el tiempo al evento independientemente de los factores clínicos?
# Comparación AIC: análoga a Tabla 6 del GLMM.
#            Pregunta: ¿mejora el PRS el ajuste del modelo de supervivencia? Criterio: ΔAIC < −2 indica mejora relevante.

# ── 4.5.1. Extraer HRs de los 5 modelos PRS ───────────────────────────────────
tabla11_datos <- lapply(names(nombres_prs), function(m) {
  extraer_prs_cox(modelos_cox_global[[m]],
                  nombre_modelo = nombres_prs[m],
                  vars_excluir  = vars_excluir_cox)
}) %>%
  do.call(rbind, .) %>%
  mutate(Modelo = factor(Modelo, levels = nombres_prs))

cat("=== HRs PRS - Cohorte FH completa (Cox) ===\n")
print(tabla11_datos %>% select(Modelo, Variable, HR, IC_95, p_valor))

# ── 4.5.2. Tabla 11: HRs del PRS ──────────────────────────────────────────────
tabla11_datos %>%
  select(Modelo, Variable, HR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Modelo") %>%
  gt_estilo(
    titulo    = "**Tabla 11. HRs del PRS - Cohorte FH completa (Cox)**",
    subtitulo = paste0("*Efecto del riesgo poligénico ajustado por el modelo base Cox · N = ",
                       cox_base$n[2], "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable PRS",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold", color = "grey30"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote  = md("HR ajustados por Edad, cLDL, cHDL, Lp(a), HTA y DM (IC 95%: exp(β ± 1,96·SE)).<br>
                    Sexo estratificado. Variables continuas estandarizadas. Frailty: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = "HR")) %>%
  guardar_gt("Tabla11_HRs_PRS_CohorteCompleta.html")

# ── 4.5.3. Tabla12. Comparación AIC - modelos PRS vs base ──────────────────────────────
modelos_cox_global_named <- setNames(modelos_cox_global, c("Base (sin PRS)", nombres_prs))
comparacion_cox_global <- comparar_modelos_cox(modelos_cox_global_named)

comparacion_cox_global %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 12. Comparación AIC - Modelos PRS Cox, Cohorte FH completa**",
    subtitulo = paste0("*AIC de la log-verosimilitud integrada · N = ", cox_base$n[2], "*")) %>%
  cols_label(
    Modelo = "Modelo",
    AIC    = "AIC",
    dAIC   = "ΔAIC") %>%
  cols_align(align = "left",  columns = Modelo) %>%
  cols_align(align = "right", columns = c(AIC, dAIC)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)) %>%
  tab_footnote(
    footnote  = md("ΔAIC = AIC_modelo − AIC_base. ΔAIC < −2 indica mejora relevante del ajuste.<br>
                    AIC calculado como −2·logLik_integrada + 2·(k + 1). C-index no reportado."),
    locations = cells_column_labels(columns = dAIC)) %>%
  guardar_gt("Tabla12_ComparacionAIC_Cox_Global.html")


# ==============================================================================
# 4.6. MODELOS PRS - ESTRATIFICADO POR SEXO (Tabla 13 + comparación AIC)
# ==============================================================================
# Replica el análisis de 4.5 por separado para mujeres y hombres.
# Sexo excluido como covariable (sin_sexo = TRUE): análisis dentro de cada sexo.

# ── 4.6.1. Ajuste de modelos por sexo ─────────────────────────────────────────
modelos_cox_m <- construir_modelos_cox(datos_surv_m, sin_sexo = TRUE)
modelos_cox_h <- construir_modelos_cox(datos_surv_h, sin_sexo = TRUE)

# ── 4.6.2. HRs del PRS por sexo ───────────────────────────────────────────────
tabla13_datos <- bind_rows(
  lapply(names(nombres_prs), function(m) {
    extraer_prs_cox(modelos_cox_m[[m]],
                    nombre_modelo = nombres_prs[m],
                    vars_excluir  = vars_excluir_cox) %>%
      mutate(Sexo = "Mujeres")
  }),
  lapply(names(nombres_prs), function(m) {
    extraer_prs_cox(modelos_cox_h[[m]],
                    nombre_modelo = nombres_prs[m],
                    vars_excluir  = vars_excluir_cox) %>%
      mutate(Sexo = "Hombres")
  })
) %>%
  mutate(Sexo  = factor(Sexo,  levels = c("Mujeres", "Hombres")),
         Modelo = factor(Modelo, levels = nombres_prs))

cat("=== HRs PRS por sexo (Cox) ===\n")
print(tabla13_datos %>% select(Sexo, Modelo, Variable, HR, IC_95, p_valor))

# ── 4.6.3. Tabla 13: HRs PRS × sexo ──────────────────────────────────────────
tabla13_datos %>%
  select(Sexo, Modelo, Variable, HR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Sexo") %>%
  gt_estilo(
    titulo    = "**Tabla 13. HRs del PRS por sexo - Cohorte FH (Cox)**",
    subtitulo = paste0("*Efecto del riesgo poligénico ajustado por el modelo base Cox, estratificado por sexo*<br>",
                       "*Mujeres N = ", modelos_cox_m$base$n[2],
                       " · Hombres N = ", modelos_cox_h$base$n[2], "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Modelo   = "Modelo PRS",
    Variable = "Variable",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Modelo, Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote  = md("HR ajustados por Edad, cLDL, cHDL, Lp(a), HTA y DM (IC 95%: exp(β ± 1,96·SE)).<br>
                    Sexo excluido como covariable. Variables continuas estandarizadas. Frailty: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = "HR")) %>%
  guardar_gt("Tabla13_HRs_PRS_Sexo.html")

# ── 4.6.4. Tabla 14. Comparación AIC por sexo ───────────────────────────────────────────
modelos_cox_m_named <- setNames(modelos_cox_m, c("Base (sin PRS)", nombres_prs))
modelos_cox_h_named <- setNames(modelos_cox_h, c("Base (sin PRS)", nombres_prs))

comparacion_cox_m <- comparar_modelos_cox(modelos_cox_m_named)
comparacion_cox_h <- comparar_modelos_cox(modelos_cox_h_named)

left_join(
  comparacion_cox_m %>% rename(AIC_M = AIC, dAIC_M = dAIC),
  comparacion_cox_h %>% rename(AIC_H = AIC, dAIC_H = dAIC),
  by = "Modelo") %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 14. Comparación AIC por sexo - Modelos PRS Cox**",
    subtitulo = paste0("*AIC de la log-verosimilitud integrada, estratificado por sexo*<br>",
                       "*Mujeres N = ", modelos_cox_m$base$n[2],
                       " · Hombres N = ", modelos_cox_h$base$n[2], "*")) %>%
  tab_spanner(
    label   = paste0("Mujeres (N = ", modelos_cox_m$base$n[2], ")"),
    columns = c(AIC_M, dAIC_M)) %>%
  tab_spanner(
    label   = paste0("Hombres (N = ", modelos_cox_h$base$n[2], ")"),
    columns = c(AIC_H, dAIC_H)) %>%
  cols_label(
    Modelo = "Modelo",
    AIC_M  = "AIC", dAIC_M = "ΔAIC",
    AIC_H  = "AIC", dAIC_H = "ΔAIC") %>%
  cols_align(align = "left",  columns = Modelo) %>%
  cols_align(align = "right", columns = c(AIC_M, dAIC_M, AIC_H, dAIC_H)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)) %>%
  tab_footnote(
    footnote  = md("ΔAIC = AIC_modelo − AIC_base por sexo. ΔAIC < −2 indica mejora relevante."),
    locations = cells_column_labels(columns = dAIC_M)) %>%
  guardar_gt("Tabla14_ComparacionAIC_Cox_Sexo.html")


# ==============================================================================
# 4.7. TEST DE INTERACCIÓN SEXO × PRS - MODELOS COX
# ==============================================================================
# Formaliza si el efecto del PRS sobre el tiempo al ECV difiere entre sexos. Análogo a la sección 3.6 del GLMM, 
# pero en escala de supervivencia.
#
# H0: el efecto del PRS sobre el hazard es igual en mujeres y hombres.
# Nota: Sexo como efecto fijo en estos modelos (no strata) para poder estimar el término de interacción. Los modelos 
# de estimación de HR (4.4-4.6) mantienen strata(Sexo). 

# Fórmula base con Sexo como efecto fijo (solo para test de interacción)
f_inter <- "Surv(t_evento, ECV_bin) ~ Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"

# ── 4.7.1. Modelos sin/con interacción para cada representación del PRS ────────
terminos_prs_inter <- c(
  grs      = "scale(GRS)",
  quintile = "Quintile",
  q5       = "Quintile_5",
  riesgo   = "Riesgo_poligenico",
  vhr      = "VHR")

modelos_cox_inter_pares <- lapply(names(terminos_prs_inter), function(nm) {
  prs <- terminos_prs_inter[[nm]]
  list(
    sin_inter = coxme(as.formula(paste(f_inter, "+", prs, "+ (1|ID_cluster)")), data = datos_surv),
    con_inter = coxme(as.formula(paste(f_inter, "+", prs, paste0("+ Sexo:", prs), "+ (1|ID_cluster)")), data = datos_surv))
})
names(modelos_cox_inter_pares) <- names(terminos_prs_inter)

# ── 4.7.2. LRT manual: loglik integrada sin vs con interacción ────────────────
lrt_cox_inter <- lapply(names(modelos_cox_inter_pares), function(nm) {
  par     <- modelos_cox_inter_pares[[nm]]
  ll_sin  <- as.numeric(par$sin_inter$loglik["Integrated"])
  ll_con  <- as.numeric(par$con_inter$loglik["Integrated"])
  chi2    <- round(-2 * (ll_sin - ll_con), 3)
  gl      <- length(par$con_inter$coefficients) - length(par$sin_inter$coefficients)
  p       <- round(pchisq(chi2, df = gl, lower.tail = FALSE), 3)
  data.frame(
    PRS   = nombres_prs[[nm]],
    Chi2  = chi2,
    gl    = gl,
    p_LRT = p,
    p_sig = p < 0.05)
}) %>% do.call(rbind, .)

cat("=== LRT Interacción Sexo × PRS (Cox) ===\n")
print(lrt_cox_inter %>% select(PRS, Chi2, gl, p_LRT))

# ── 4.7.3. Tabla 15: Test de interacción Cox ──────────────────────────────────
lrt_cox_inter %>%
  select(PRS, Chi2, gl, p_LRT, p_sig) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 15. Test de interacción Sexo × PRS - Modelos Cox**",
    subtitulo = paste0("*Likelihood Ratio Test · N = ", cox_base$n[2], " · Sexo como efecto fijo en modelos de interacción*")) %>%
  cols_label(
    PRS   = "Representación PRS",
    Chi2  = "Chi²",
    gl    = "gl",
    p_LRT = "Valor p") %>%
  cols_hide("p_sig") %>%
  cols_align(align = "left",  columns = PRS) %>%
  cols_align(align = "right", columns = c(Chi2, gl, p_LRT)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p_LRT, rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("LRT: Chi² = -2·(logLik_sin - logLik_con) ~ χ² con gl = nº términos de interacción.<br>
                    Modelos de interacción: Sexo incluido como efecto fijo (no strata) para permitir la estimación<br>
                    del término Sexo × PRS. Frailty: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = Chi2)) %>%
  tab_source_note(
    source_note = md(paste0(
      "Cohorte FH: N = ", cox_base$n[2], " · Ajustado por Sexo, Edad, cLDL, cHDL, Lp(a), HTA, DM + (1|ID_cluster)"))) %>%
  guardar_gt("Tabla15_Interaccion_SexoPRS_Cox.html")

# ==============================================================================
# 4.8. ANÁLISIS DE SENSIBILIDAD - EVENTOS INCIDENTES (Tabla 16)
# ==============================================================================
# Objetivo: verificar que los resultados del modelo Cox no están distorsionados por el sesgo de left truncation derivado
# de los eventos prevalentes (ECV ocurrido antes de la inclusión en la cohorte).
#
# Estrategia: restringir el análisis a eventos incidentes (post-inclusión) mediante t_post = Edad_ECV − Edad_inclusion 
# (ya calculado en 4.1). Individuos con ECV previo a la inclusión pasan a ser censurados en t=0+.
# Si el GRS sigue siendo NS en este subconjunto, el resultado principal no está contaminado por sesgo de supervivencia selectiva.

# ── 4.8.1. Dataset de eventos incidentes ──────────────────────────────────────
# Excluir individuos con ECV previa a la inclusión (evento prevalente): su t_post sería <= 0, no aportan información de incidencia
datos_surv_inc <- datos_surv %>%
  filter(!(ECV_bin == 1 & t_post <= 0)) %>%
  mutate(ECV_inc = case_when(
    ECV_bin == 1 & t_post > 0 ~ 1L,   # evento incidente confirmado
    TRUE                      ~ 0L))  # censurado (sin evento o evento prevalente)

cat("=== Dataset sensibilidad - eventos incidentes ===\n")
cat("N total:                  ", nrow(datos_surv_inc),          "\n")
cat("Eventos incidentes:       ", sum(datos_surv_inc$ECV_inc),   "\n")
cat("Censurados (incl. prev.): ", sum(datos_surv_inc$ECV_inc==0),"\n")
cat("Excluidos (t_post <= 0):  ", nrow(datos_surv) - nrow(datos_surv_inc), "\n")

# ── 4.8.2. Modelo Cox GRS - eventos incidentes ────────────────────────────────
cox_inc_grs <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                       scale(LpA_0) + HTA_bin + DM_bin + scale(GRS) + (1 | ID_cluster),
                     data = datos_surv_inc)
summary(cox_inc_grs)

# ── 4.8.3. Extraer HR del GRS ─────────────────────────────────────────────────
hr_grs_inc <- extraer_hrs(cox_inc_grs)["scale(GRS)", ]

cat("\n=== HR del GRS — eventos incidentes ===\n")
cat("HR  =", hr_grs_inc$HR, "\n")
cat("IC95= (", hr_grs_inc$IC_2.5, "–", hr_grs_inc$IC_97.5, ")\n")
cat("p   =", hr_grs_inc$p_valor, "\n")

# ── 4.8.4. Tabla 16: Sensibilidad ─────────────────────────────────────────────
extraer_hrs(cox_inc_grs) %>%
  tibble::rownames_to_column("Var_raw") %>%
  mutate(
    Variable = ifelse(Var_raw %in% names(var_labels_base_cox), var_labels_base_cox[Var_raw],
                      ifelse(Var_raw == "scale(GRS)", "GRS continuo (por DE)", Var_raw)),
    IC_95    = paste0("(", IC_2.5, " – ", IC_97.5, ")"),
    p_sig    = p_valor == "<0.001" |
      (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05),
    PRS      = Var_raw == "scale(GRS)") %>%
  select(Variable, HR, IC_95, p_valor, p_sig, PRS) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 16. Análisis de sensibilidad - Eventos incidentes (Cox)**",
    subtitulo = paste0("*Restricción a ECV post-inclusión · N = ", nrow(datos_surv_inc),
                       " · Eventos incidentes = ", sum(datos_surv_inc$ECV_inc), "*")) %>%
  cols_hide(c("p_sig", "PRS")) %>%
  cols_label(
    Variable = "Variable",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("HR ajustados (IC 95%: exp(β ± 1,96·SE)). Variables continuas estandarizadas.<br>
                    Tiempo al evento: t_post = Edad_ECV − Edad_inclusión (desde la inclusión).<br>
                    Individuos con ECV previa a la inclusión reclasificados como censurados."),
    locations = cells_column_labels(columns = "HR")) %>%
  tab_source_note(
    source_note = md(paste0(
      "Comparar con modelo principal (Tabla 10): GRS HR = ", extraer_hrs(cox_base)["scale(GRS)", "HR"],
      " (p = ", extraer_hrs(cox_base)["scale(GRS)", "p_valor"], ")."))) %>%
  guardar_gt("Tabla16a_Sensibilidad_EventosIncidentes.html")

# ── 4.8.2b. Modelos adicionales - eventos incidentes ──────────────────────────
cox_inc_q5     <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                          scale(LpA_0) + HTA_bin + DM_bin + Quintile_5 + (1|ID_cluster), data = datos_surv_inc)

cox_inc_quintile     <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                                scale(LpA_0) + HTA_bin + DM_bin + Quintile + (1|ID_cluster), data = datos_surv_inc)

cox_inc_riesgo <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                          scale(LpA_0) + HTA_bin + DM_bin + Riesgo_poligenico + (1|ID_cluster), data = datos_surv_inc)

# Estratificado por hombres 
datos_surv_inc_h <- filter(datos_surv_inc, Sexo == "Hombre")

cox_inc_q5_h   <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                          HTA_bin + DM_bin + Quintile_5 + (1|ID_cluster), data = datos_surv_inc_h)

cox_inc_quintile_h   <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                                HTA_bin + DM_bin + Quintile + (1|ID_cluster), data = datos_surv_inc_h)

cox_inc_riesgo_h <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                            HTA_bin + DM_bin + Riesgo_poligenico + (1|ID_cluster), data = datos_surv_inc_h)

cat("=== Eventos incidentes — HRs PRS adicionales ===\n")
cat("Q5 binario (global):  ");  print(extraer_hrs(cox_inc_q5)["Quintile_5Si", ])
cat("Quintiles (global):  ");  print(extraer_hrs(cox_inc_quintile)[c("QuintileQ1", "QuintileQ2", "QuintileQ3", "QuintileQ4", "QuintileQ5"), ])
cat("Riesgo poligénico (global):\n"); print(extraer_hrs(cox_inc_riesgo)[c("Riesgo_poligenicoIntermedio","Riesgo_poligenicoAlto"), ])
cat("Q5 binario (hombres): "); print(extraer_hrs(cox_inc_q5_h)["Quintile_5Si", ])
cat("Quintiles (hombres):  ");  print(extraer_hrs(cox_inc_quintile_h)[c("QuintileQ1", "QuintileQ2", "QuintileQ3", "QuintileQ4", "QuintileQ5"), ])
cat("Riesgo poligénico (hombres):\n"); print(extraer_hrs(cox_inc_riesgo_h)[c("Riesgo_poligenicoIntermedio","Riesgo_poligenicoAlto"), ])

# ── 4.8.3b. Tabla 16 ampliada: GRS + PRS categórico global y hombres ──────────
bind_rows(
  # GRS global
  extraer_hrs(cox_inc_grs) %>% rownames_to_column("Var_raw") %>% filter(Var_raw == "scale(GRS)") %>%
    mutate(Grupo = "Cohorte global", Variable = "GRS continuo (por DE)"),
  # Q5 global
  extraer_hrs(cox_inc_q5) %>% rownames_to_column("Var_raw") %>% filter(Var_raw == "Quintile_5Si") %>%
    mutate(Grupo = "Cohorte global", Variable = "Quintil 5 vs. Q1-Q4"),
  # Riesgo global
  extraer_hrs(cox_inc_riesgo) %>% rownames_to_column("Var_raw") %>% filter(grepl("Riesgo_poligenico", Var_raw)) %>%
    mutate(Grupo    = "Cohorte global",
           Variable = recode(Var_raw, "Riesgo_poligenicoIntermedio" = "Riesgo Intermedio vs. Bajo",
                             "Riesgo_poligenicoAlto"       = "Riesgo Alto vs. Bajo")),
  # Q5 hombres
  extraer_hrs(cox_inc_q5_h) %>% rownames_to_column("Var_raw") %>% filter(Var_raw == "Quintile_5Si") %>%
    mutate(Grupo = "Hombres FH", Variable = "Quintil 5 vs. Q1-Q4"),
  # Riesgo hombres
  extraer_hrs(cox_inc_riesgo_h) %>% rownames_to_column("Var_raw") %>% filter(grepl("Riesgo_poligenico", Var_raw)) %>%
    mutate(Grupo    = "Hombres FH",
           Variable = recode(Var_raw, "Riesgo_poligenicoIntermedio" = "Riesgo Intermedio vs. Bajo",
                             "Riesgo_poligenicoAlto"       = "Riesgo Alto vs. Bajo"))) %>%
  mutate(
    IC_95 = paste0("(", IC_2.5, " - ", IC_97.5, ")"),
    p_sig = p_valor == "<0.001" |
      (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)) %>%
  select(Grupo, Variable, HR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Grupo") %>%
  gt_estilo(
    titulo    = "**Tabla 16b. Análisis de sensibilidad - Eventos incidentes (Cox)**",
    subtitulo = paste0("*Restricción a ECV post-inclusión · N = ", nrow(datos_surv_inc),
                       " · Eventos incidentes = ", sum(datos_surv_inc$ECV_inc), "*")) %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable PRS",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p") %>%
  cols_align(align = "left",  columns = c(Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote  = md(paste0("HR ajustados (IC 95%: exp(β ± 1,96·SE)). Variables continuas estandarizadas.<br>
                    t_post = Edad_ECV − Edad_inclusión. Individuos con ECV previa reclasificados como censurados.<br>
                    Hombres: n = ", nrow(datos_surv_inc_h), ". Eventos incidentes = ", sum(datos_surv_inc_h$ECV_inc))),
    locations = cells_column_labels(columns = "HR")) %>%
  guardar_gt("Tabla16b_Sensibilidad_EventosIncidentes_Ampliada.html")








