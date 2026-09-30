######################################################################################################################
# TFM: Impacto del Riesgo Poligénico Cardiovascular en Hipercolesterolemia Familiar
# Autora:     Tamara Noya Mosquera
# Directores: Dra. Teresa Padró - Institut de Recerca Sant Pau (IR Sant Pau)
#             Ariel Ernesto Cariaga Martínez - Universidad Internacional de la Rioja (UNIR)
# Programa:   Máster en Bioinformática, UNIR (curso 2025-2026)
#
# Dependencias: readxl, dplyr, tidyr, anytime, gtsummary, gt, ggplot2, ggsignif, lme4, performance, pROC, tibble, 
#               survival, survminer, coxme, scales, patchwork
######################################################################################################################


######################################################################################################################
# BLOQUE 1 - CONFIGURACIÓN, CARGA Y PREPARACIÓN DE DATOS
######################################################################################################################

# ==============================================================================
# 1.1. ENTORNO DE TRABAJO
# ==============================================================================

setwd("C:/Users/LG/Documents/BIOINFO UNIR/TFM/Analisis_R/Datos")

# Directorio de salida para tablas y figuras 
ruta_out <- "C:/Users/LG/Documents/BIOINFO UNIR/TFM/Analisis_R/Figuras/Figuras_R_TMF"

# Semilla de reproducibilidad - garantiza resultados idénticos en cada ejecución.
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
library(lme4)        # Modelos lineales/generalizados mixtos (glmer - GLMM)
library(performance) # Diagnóstico de modelos mixtos (VIF, colinealidad)
library(pROC)        # Curvas ROC, AUC y comparación de modelos (DeLong)
library(tibble)      # rownames_to_column y manipulación de tibbles
library(survival)    # Análisis de supervivencia: Surv(), coxph(), survfit()
library(survminer)   # Visualización de curvas Kaplan-Meier (ggsurvplot)
library(coxme)       # Modelos de Cox con frailty log-normal (estructura familiar)
library(scales)      # Formatos de ejes en ggplot 
library(patchwork)   # Composición de múltiples gráficos ggplot en una figura
library(car)         # VIF (Factor de Inflación de Varianza): diagnóstico de colinealidad

# Tema compacto para tablas gtsummary (reduce el espaciado por defecto)
theme_gtsummary_compact()

# ==============================================================================
# 1.3. CONSTANTES DE ESTILO
# ==============================================================================

# ── Paletas de colores ─────────────────────────────────────────────────────────
# Grupo principal (FH vs No FH)
colores_fh    <- c("No FH" = "#B2BABB", "FH" = "#1F5299")
# Grupo de estudio (3 grupos: No FH / FH Sin Evento / FH Con Evento)
colores_grupo <- c("No FH" = "#B2BABB", "FH No Evento" = "#8FA4B5", "FH Evento" = "#3182CE")
# Subgrupos FH por evento (para gráficos dentro de la cohorte FH)
colores_ecv   <- c("FH No Evento" = "#8FA4B5", "FH Evento" = "#3182CE")
# Fenotipo de ECV (Sin Evento / Precoz <65a / Tardío ≥65a)
colores_fenotipo <- c("Sin Evento" = "#8FA4B5", "Precoz" = "#5DADE2", "Tardío" = "#2B6CB0")
# Quintiles del PRS (gradiente de riesgo bajo → alto)
colores_quintile <- c("Q1" = "#4CB97B", "Q2" = "#FAD390", "Q3" = "#F89F73", "Q4" = "#D35400", "Q5" = "#A83232")
# Riesgo poligénico agrupado (Q1=Bajo, Q2-Q4=Intermedio, Q5=Alto)
colores_riesgo <- c("Bajo" = "#4CB97B", "Intermedio" = "#F89F73", "Alto" = "#A83232")
# Q5 vs resto (análisis binario de riesgo alto)
colores_q5  <- c("Q1-Q4" = "#80CFA3", "Q5" = "#A83232")
# VHR (Very High Risk: subconjunto de Q5 con criterios adicionales)
colores_vhr <- c("No" = "#80CFA3", "Si" = "#8B1E1E")
# Genotipos de SNPs (dosage 0/1/2)
colores_snp <- c("Hom. referencia (0)" = "#ECEFF1", "Heterocigoto (1)" = "#90A4AE", "Hom. riesgo (2)" = "#455A64")

# ── Variables y funciones auxiliares para figuras ───────────────────────────────────────────────────────────
# Aplica un tema bas que se combina con theme() adicional en cada gráfico para ajustes específicos.
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
# BLOQUE 2 - ANÁLISIS DESCRIPTIVO
######################################################################################################################
# Objetivo: caracterizar la cohorte y comparar grupos (FH vs No FH; Sin/Con Evento)

# ==============================================================================
# 2.0. CONSTANTES COMPARTIDAS DEL BLOQUE
# ==============================================================================
# Centralizar aquí variables, etiquetas y estadísticos reutilizados en varias tablas

# Variables analíticas principales 
vars_t1 <- c("Edad_inclusion", "Sexo", "CT_0", "cLDL_0", "cHDL_0", "TG_0", "ApoB_0", "LpA_0", "HTA_bin", "DM_bin",
             "TipAleloRLDL_0", "MutApoB_0", "LLT_intensidad_0", "AñosTtoEstatinas_0",
             "GRS", "Quintile", "VHR")

# Etiquetas legibles para todas las tablas 
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
# Tabla principal del análisis - compara los dos grupos FH según el outcome primario (evento cardiovascular).
# Test: Wilcoxon rank-sum (continuas), Fisher exact (categóricas).
# Justificación Fisher: algunas celdas tienen n<5 (MutApoB_0, VHR).
tabla1 <- datos_fh %>%
  tbl_summary(
    by        = Grupo_estudio,
    include   = all_of(vars_t1),
    statistic = estadisticos_base,
    digits    = list(all_continuous() ~ 1, GRS ~ 3),
    missing   = "no",     
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

# ── Parte B: subgrupos FH - descriptivo únicamente, sin p ─────────────────────
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
      missing   = "no",            
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
# Tres figuras complementarias que documentan la  asociación del GRS con el fenotipo cardiovascular en 
# la cohorte FH y en comparación con No FH.
#
# Estructura de análisis: Global → pairwise 3 grupos → por sexo (global 3 grupos + pairwise)
#
# Diseño:
#    - Gráfico A: visión binaria No FH vs FH 
#    - Gráfico B: 3 grupos con todos los pares 
#    - Gráfico C: estratificado por sexo (solo FH (n No FH por sexo demasiado pequeña))
#
# Tests: Wilcoxon rank-sum con corrección de Bonferroni (3 comparaciones en B)
# Se muestran todas las comparaciones incluyendo "ns"

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
# Evalúa si el PRS discriminala precocidad de ECV.
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
cat("Global — KW:               p =", round(kw_grs_feno$p.value,   3), "\n")
cat("Global — Sin Ev vs Precoz: p =", round(p_sine_precoz,         3), "(Bonferroni)\n")
cat("Global — Sin Ev vs Tardío: p =", round(p_sine_tardio,         3), "(Bonferroni)\n")
cat("Global — Precoz vs Tardío: p =", round(p_prec_tardio,         3), "(Bonferroni)\n")
cat("Mujeres — KW:              p =", round(kw_grs_feno_m$p.value, 3), "\n")
cat("Mujeres — Sin vs Precoz:   p =", round(pw_grs_feno_m["Precoz","Sin Evento"], 3), "(Bonferroni)\n")
cat("Mujeres — Sin vs Tardío:   p =", round(pw_grs_feno_m["Tardío","Sin Evento"], 3), "(Bonferroni)\n")
cat("Mujeres — Precoz vs Tardío:p =", round(pw_grs_feno_m["Tardío","Precoz"],    3), "(Bonferroni)\n")
cat("Hombres — KW:              p =", round(kw_grs_feno_h$p.value, 3), "\n")
cat("Hombres — Sin vs Precoz:   p =", round(pw_grs_feno_h["Precoz","Sin Evento"], 3), "(Bonferroni)\n")
cat("Hombres — Sin vs Tardío:   p =", round(pw_grs_feno_h["Tardío","Sin Evento"], 3), "(Bonferroni)\n")
cat("Hombres — Precoz vs Tardío:p =", round(pw_grs_feno_h["Tardío","Precoz"],    3), "(Bonferroni)\n")

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
#   · FH Sin vs FH Con Evento: ¿PRS predice ECV dentro de FH? 
#
# Función chi_analisis(): ejecuta los 4 tests (global + 3 pairwise) para cualquier variable y dataset, evitando repetición de código.
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
# Subtítulo: p global de 3 grupos por sexo (incl. No FH) + p FH Sin vs FH Con.
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
# El chi² FH Sin vs FH Con en hombres (p=0.002) es estadísticamente significativo, pero no identifica qué quintil(es) 
# son responsables. Este análisis post-hoc descompone el efecto para interpretar el patrón de asociación.
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
# Pregunta B: ¿tiene cada quintil más o menos eventos de los esperados dado el 63.1% global? 
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
    ci_high  = round(binom.test(n_cone, n_total)$conf.int[2] * 100, 1)) %>%
  ungroup()

cat("\n=== Tasa de eventos por quintil vs esperada (", round(p_eventos_h*100,1), "%) ===\n")
print(binom_quintile_h %>% select(Quintile, n_total, pct_cone, ci_low, ci_high))

# ── 2.7.4. Gráfico: tasa de eventos por quintil con IC 95% ─────────────────
# Visualiza el patrón no monotónico. La línea discontinua marca la tasa global esperada (63.1%).
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
  geom_text(aes(label = paste0(pct_cone, "%"),
                y     = ci_high + 2),
            size = 2.7, color = "grey20") +
  scale_fill_manual(values = colores_quintile) +
  scale_y_continuous(limits = c(0, 90), labels = function(x) paste0(x, "%"), breaks = seq(0, 80, 20)) +
  labs(
    title    = "Tasa de evento cardiovascular por quintil - Hombres FH",
    subtitle = paste0("Línea discontinua = tasa media en hombres FH (", round(p_eventos_h * 100, 1), "%).  Barras = IC 95% binomial\n",
                      "Fisher pairwise vs. Q1 (Bonferroni): Q2 p=", pairwise_q_h$p_bonferroni[1],
                      "; Q3 p=", pairwise_q_h$p_bonferroni[2],
                      "; Q4 p=", pairwise_q_h$p_bonferroni[3],
                      "; Q5 p=", pairwise_q_h$p_bonferroni[4]),
    caption  = paste0("Análisis post-hoc exploratorio · Hombres FH N=", nrow(datos_fh_h),
                      ". Fisher exacto con corrección de Bonferroni (4 comparaciones: Q2-Q5 vs. Q1)"),
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

# Variables del modelo base
var_labels_base <- c(
  "SexoHombre"             = "Sexo (Hombre vs. Mujer)",
  "scale(Edad_inclusion)"  = "Edad en la inclusión (por DE)",
  "scale(cLDL_0)"          = "cLDL basal (por DE)",
  "scale(cHDL_0)"          = "cHDL basal (por DE)",
  "scale(LpA_0)"           = "Lp(a) basal (por DE)",
  "HTA_binSi"              = "Hipertensión arterial (Sí vs. No)",
  "DM_binSi"               = "Diabetes mellitus (Sí vs. No)")

# Variables PRS 
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
# 3.1. ANÁLISIS DEL EFECTO FAMILIAR - ICC Y DESIGN EFFECT
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
    titulo    = "**Tabla ICC. Análisis de estructura familiar - Justificación del GLMM**") %>%
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
# 3.2. ANÁLISIS UNIVARIANTE - OR NO AJUSTADOS (Tabla 3)
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
# 3.3. MODELO BASE GLMM - COHORTE FH COMPLETA (Tabla 4)
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
# Se guardan en objetos para reutilizar 
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
# 3.4. MODELOS PRS - COHORTE FH COMPLETA 
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
# 3.5. MODELOS PRS - ESTRATIFICADO POR SEXO 
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
                       "*Mujeres N = ", nobs(modelos_m$base), " · Hombres N = ", nobs(modelos_h$base), "*")) %>%
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
                       "*Mujeres N = ", nobs(modelos_m$base), " · Hombres N = ", nobs(modelos_h$base), "*")) %>%
  tab_spanner(label   = paste0("Mujeres (N = ", nobs(modelos_m$base), ")"),
              columns = c(AIC_M, AUC_M, dAIC_M, dAUC_M)) %>%
  tab_spanner(label   = paste0("Hombres (N = ", nobs(modelos_h$base), ")"),
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
# BLOQUE 4 - ANÁLISIS DE SUPERVIVENCIA (KAPLAN-MEIER + COX CON FRAILTY)
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
# determinar su tiempo de censura. Para el resto: tiempo al evento = Edad_ECV - 18 (si ECV) o Edad_obs - 18 (si censurado)
# Se restan 18 años porque la HF es una enfermedad congénita con riesgo cardiovascular acumulado desde la infancia; 
# anclar en los 18 años captura la exposición de riesgo desde la edad adulta temprana 
#
# t_post (desde la inclusión): reservado para el análisis de sensibilidad.
# Permite restringir el análisis a eventos incidentes (post-inclusión), que descarta el sesgo de left truncation 
# por eventos prevalentes.

datos_surv <- datos_fh %>%
  filter(!(Muerte == "Si" & is.na(Edad_muerte) & ECV != "Si")) %>%
  mutate(
    # Tiempo principal: desde los 18 años
    t_evento = case_when(
      ECV_bin == 1 ~ Edad_ECV       - 18,
      TRUE         ~ Edad_obs       - 18),
    # Tiempo desde inclusión: para análisis de sensibilidad (4.8)
    t_post   = case_when(
      ECV_bin == 1 ~ Edad_ECV       - Edad_inclusion,
      TRUE         ~ Edad_obs       - Edad_inclusion))

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
# Los modelos finales incluyen el frailty term (1 | ID_cluster).

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
# Figura 2B: KM × Quintile en hombres (5 grupos). Pregunta: ¿se observa el patrón no monotónico (Q2 > Q1 ≈ Q5)?
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
# Visualiza el patrón no monotónico en hombres: si Q2 tiene peor supervivencia que el resto y Q5 es similar a Q1
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
# 4.4. MODELO COX BASE CON FRAILTY 
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
# La discriminación del modelo GLMM (AUC = 0.870) proporciona la medida de discriminación principal. Los modelos Cox 
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
# 4.5. MODELOS PRS - COHORTE FH COMPLETA 
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
# 4.6. MODELOS PRS - ESTRATIFICADO POR SEXO 
# ==============================================================================
# Replica el análisis de 4.5 por separado para mujeres y hombres.
# El GLMM mostró dimorfismo sexual significativo. Se evalúa si el patrón no monotónico en hombres persiste en la escala temporal.
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
# 4.8. ANÁLISIS DE SENSIBILIDAD - EVENTOS INCIDENTES 
# ==============================================================================
# Objetivo: verificar que los resultados del modelo Cox no están distorsionados por el sesgo de left truncation derivado
# de los eventos prevalentes (ECV ocurrido antes de la inclusión).
#
# Estrategia: restringir el análisis a eventos incidentes (post-inclusión) filtrando sujetos con t_post > 0 o sin evento.
# Los individuos con ECV previo a la inclusión (evento prevalente) quedan excluidos de la muestra analítica.

# ── 4.8.1. Dataset de eventos incidentes ──────────────────────────────────────
# Mantener individuos sanos (ECV_bin == 0) y aquellos con eventos incidentes (t_post > 0)
datos_surv_inc <- datos_surv %>%
  filter(t_post > 0 | ECV_bin == 0) %>%
  mutate(ECV_inc = if_else(ECV_bin == 1, 1L, 0L))

# Subset específico para el estrato de hombres
datos_surv_inc_h <- datos_surv_inc %>% 
  filter(Sexo == "Hombre")

cat("=== Dataset sensibilidad - eventos incidentes ===\n")
cat("N total global:      ", nrow(datos_surv_inc), "\n")
cat("Eventos incidentes:  ", sum(datos_surv_inc$ECV_inc), "\n")
cat("N total hombres:     ", nrow(datos_surv_inc_h), "\n")
cat("Eventos inc. hombres:", sum(datos_surv_inc_h$ECV_inc), "\n")

# ── 4.8.2. Modelos Cox GRS - eventos incidentes ────────────────────────────────
cox_inc_grs <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                       scale(LpA_0) + HTA_bin + DM_bin + scale(GRS) + (1 | ID_cluster), data = datos_surv_inc)

cox_inc_q5     <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                          scale(LpA_0) + HTA_bin + DM_bin + Quintile_5 + (1 | ID_cluster), data = datos_surv_inc)

cox_inc_riesgo <- coxme(Surv(t_post, ECV_inc) ~ strata(Sexo) + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + 
                          scale(LpA_0) + HTA_bin + DM_bin + Riesgo_poligenico + (1 | ID_cluster), data = datos_surv_inc)

# Modelos en el estrato de hombres
cox_inc_grs_h  <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                          HTA_bin + DM_bin + scale(GRS) + (1 | ID_cluster), data = datos_surv_inc_h)

cox_inc_q5_h   <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                          HTA_bin + DM_bin + Quintile_5 + (1 | ID_cluster), data = datos_surv_inc_h)

cox_inc_quintile_h <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                              HTA_bin + DM_bin + Quintile + (1 | ID_cluster), data = datos_surv_inc_h)

cox_inc_riesgo_h <- coxme(Surv(t_post, ECV_inc) ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +
                            HTA_bin + DM_bin + Riesgo_poligenico + (1 | ID_cluster), data = datos_surv_inc_h)

# ── 4.8.3. Estructuración de filas para Tabla 16 ────────────────────
grupo_global <- paste0("Cohorte global (N = ", cox_inc_grs$n[2],   ". Eventos = ", cox_inc_grs$n[1],   ")")
grupo_h      <- paste0("Hombres FH (n = ",     cox_inc_grs_h$n[2], ". Eventos = ", cox_inc_grs_h$n[1], ")")

# Función interna auxiliar para homogeneizar las extracciones del PRS
extraer_prs_inc <- function(modelo, patron_var, etiqueta, grupo) {
  extraer_hrs(modelo) %>%
    rownames_to_column("Var_raw") %>%
    filter(grepl(patron_var, Var_raw)) %>%
    mutate(Grupo = grupo, Variable = etiqueta)
}

# ── 4.8.4. Tabla 16 ──────────────────────────────────
tabla16_data <- bind_rows(
  # Variables clínicas (Cohorte global)
  extraer_hrs(cox_inc_grs) %>%
    rownames_to_column("Var_raw") %>%
    filter(Var_raw %in% names(var_labels_base_cox)) %>%
    mutate(Grupo = grupo_global, Variable = var_labels_base_cox[Var_raw]),
  
  # 2. Representaciones del PRS (Cohorte global)
  extraer_prs_inc(cox_inc_grs,    "scale\\(GRS\\)", "GRS continuo (por DE)",      grupo_global),
  extraer_prs_inc(cox_inc_q5,     "Quintile_5Si",   "Quintil 5 vs. Q1–Q4",        grupo_global),
  extraer_prs_inc(cox_inc_riesgo, "Intermedio",     "Riesgo Intermedio vs. Bajo", grupo_global),
  extraer_prs_inc(cox_inc_riesgo, "Alto",           "Riesgo Alto vs. Bajo",       grupo_global),
  
  # 3. Representaciones del PRS (Estrato Hombres FH)
  extraer_prs_inc(cox_inc_grs_h,      "scale\\(GRS\\)", "GRS continuo (por DE)",      grupo_h),
  extraer_prs_inc(cox_inc_quintile_h, "QuintileQ2",     "Quintil 2 vs. Q1",           grupo_h),
  extraer_prs_inc(cox_inc_quintile_h, "QuintileQ3",     "Quintil 3 vs. Q1",           grupo_h),
  extraer_prs_inc(cox_inc_quintile_h, "QuintileQ4",     "Quintil 4 vs. Q1",           grupo_h),
  extraer_prs_inc(cox_inc_quintile_h, "QuintileQ5",     "Quintil 5 vs. Q1",           grupo_h),
  extraer_prs_inc(cox_inc_q5_h,       "Quintile_5Si",   "Quintil 5 vs. Q1–Q4",        grupo_h),
  extraer_prs_inc(cox_inc_riesgo_h,   "Intermedio",     "Riesgo Intermedio vs. Bajo", grupo_h),
  extraer_prs_inc(cox_inc_riesgo_h,   "Alto",           "Riesgo Alto vs. Bajo",       grupo_h)
)

# Renderizado y exportación
tabla16_gt <- tabla16_data %>%
  mutate(
    IC_95 = paste0("(", IC_2.5, " – ", IC_97.5, ")"),
    p_sig = p_valor == "<0.001" |
      (!is.na(suppressWarnings(as.numeric(p_valor))) & suppressWarnings(as.numeric(p_valor)) < 0.05)
  ) %>%
  select(Grupo, Variable, HR, IC_95, p_valor, p_sig) %>%
  gt(groupname_col = "Grupo") %>%
  gt_estilo(
    titulo    = "**Tabla 16. Análisis de sensibilidad - Eventos incidentes (Cox)**",
    subtitulo = "*Restricción a ECV post-inclusión*") %>%
  cols_hide("p_sig") %>%
  cols_label(
    Variable = "Variable",
    HR       = "HR",
    IC_95    = "IC 95%",
    p_valor  = "Valor p"
  ) %>%
  cols_align(align = "left",  columns = c(Variable, IC_95)) %>%
  cols_align(align = "right", columns = c(HR, p_valor)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_valor", rows = p_sig)
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_footnote(
    footnote  = md("HR ajustados mediante modelos de Cox con frailty familiar (IC 95%: exp(β ± 1,96·SE)). Variables continuas estandarizadas.<br>
                    Sujetos con evento cardiovascular previo a la inclusión fueron excluidos del análisis de sensibilidad.<br>
                    Modelos globales ajustados por strata(Sexo) + Edad, cLDL, cHDL, Lp(a), HTA, DM + (1|ID_cluster).<br>
                    Modelos de hombres estiman el efecto dentro del estrato (específico por sexo). t_post = Edad_ECV - Edad_inclusión."),
    locations = cells_column_labels(columns = "HR")
  )
guardar_gt(tabla16_gt, "Tabla16_Sensibilidad_EventosIncidentes.html")


######################################################################################################################
# BLOQUE 5 - ANÁLISIS DE SNPs INDIVIDUALES
######################################################################################################################
#
# Objetivo: evaluar el efecto de cada uno de los 12 SNPs del score de forma independiente, sin asumir los pesos de la población general
#
# Justificación:
#   - Los pesos del score fueron calibrados en población general, no en FH monogénica
#   - Algunos SNPs pueden tener efectos diferenciales en contexto LDLR+, especialmente considerando dimorfismo 
#     sexual 
#   - El análisis individual permite identificar qué variantes, si alguna, se asocian con ECV en FH independientemente 
#     del score compuesto


# ==============================================================================
# 5.0. CONSTANTES Y FUNCIONES AUXILIARES DEL BLOQUE
# ==============================================================================
# ── Definición de SNPs, genes y alelos de riesgo ──────────────────────────────
snp_raw <- c("rs10507391", "rs17222842", "rs9315050",  "rs17216473", "rs1333049",  "rs9982601",  "rs10455872_LpA",
             "rs17464857", "rs9818870",  "rs12526453", "rs6725887",  "rs501120")

snp_genes <- c("ALOX5AP", "ALOX5AP", "ALOX5AP", "ALOX5AP", "CDKN2A/B", "SLC5A3/KCNE2", "LPA", "TAF1A",
               "MRAS", "PHACTR1", "WDR12", "CXCL12")

snp_risk_alleles <- c(
  rs10507391     = "A",  rs17222842     = "G",  rs9315050      = "A",
  rs17216473     = "A",  rs1333049      = "C",  rs9982601      = "T",
  rs10455872_LpA = "G",  rs17464857     = "T",  rs9818870      = "T",
  rs12526453     = "C",  rs6725887      = "C",  rs501120       = "T")

snp_nombres <- gsub("_LpA", "", snp_raw)   # nombres limpios para tablas y figuras

# ── normalizar_geno(): estandariza la codificación de genotipos ───────────────
# Ordena alfabéticamente los dos alelos (AG = GA → "AG") para garantizar representación única independientemente 
# del orden de lectura del genotipo.
normalizar_geno <- function(x) {
  factor(sapply(as.character(x), function(g) {
    if (is.na(g)) return(NA_character_)
    paste(sort(strsplit(g, "")[[1]]), collapse = "")
  }))
}

# ── crear_dosage(): convierte genotipo a número de alelos de riesgo (0/1/2) ───
# 0 = homocigoto referencia, 1 = heterocigoto, 2 = homocigoto riesgo.
crear_dosage <- function(snp_var, alelo_riesgo, df) {
  geno <- as.character(df[[snp_var]])
  sapply(strsplit(geno, ""), function(a) {
    if (any(is.na(a))) return(NA_integer_)
    as.integer(sum(a == alelo_riesgo))
  })
}

# ── calcular_freq_snp() ───────────────────────────────────────────────────────
# Calcula frecuencias genotípicas (Hom.ref / Het / Hom.riesgo) y RAF para los 12 SNPs estratificado por cualquier 
# variable de agrupación del dataset. Incluye chi-cuadrado con simulación Monte Carlo entre los niveles solicitados.
#
# Uso:
#   calcular_freq_snp(datos,    "Grupo_estudio", c("No FH","FH No Evento","FH Evento"))
#   calcular_freq_snp(datos_fh, "Sexo",          c("Mujer","Hombre"))
#   calcular_freq_snp(datos_fh, "Fenotipo_ECV",  c("Sin Evento","Precoz","Tardío"))
calcular_freq_snp <- function(df, var_grupo, niveles_grupo) {
  do.call(rbind, lapply(seq_along(snp_raw), function(i) {
    dose_var <- dose_vars[i]
    filas <- do.call(rbind, lapply(niveles_grupo, function(g) {
      mask <- !is.na(df[[var_grupo]]) & as.character(df[[var_grupo]]) == g
      dose <- df[[dose_var]][mask]
      n    <- sum(!is.na(dose))
      if (n == 0) return(NULL)
      data.frame(
        Grupo        = g,
        SNP          = snp_nombres[i],
        Gen          = snp_genes[i],
        N            = n,
        Hom_ref_pct  = round(sum(dose == 0, na.rm=TRUE) / n * 100, 1),
        Het_pct      = round(sum(dose == 1, na.rm=TRUE) / n * 100, 1),
        Hom_risk_pct = round(sum(dose == 2, na.rm=TRUE) / n * 100, 1),
        RAF          = round((2*sum(dose==2,na.rm=TRUE) + sum(dose==1,na.rm=TRUE)) / (2*n), 4))
    }))
    cols_ok  <- niveles_grupo[niveles_grupo %in% colnames(table(as.character(df[[snp_raw[i]]]), as.character(df[[var_grupo]])))]
    tab <- table(as.character(df[[snp_raw[i]]]), as.character(df[[var_grupo]]))[, cols_ok, drop = FALSE]
    filas$p_chi <- tryCatch(
      round(chisq.test(tab, simulate.p.value=TRUE, B=2000)$p.value, 4),
      error = function(e) NA_real_)
    filas
  })) %>% {rownames(.) <- NULL; .}
}
# ── graficar_freq_snp() ───────────────────────────────────────────────────────
# Gráfico de barras apiladas estándar a partir de una tabla de calcular_freq_snp().
# etiquetas_x: named vector para renombrar el eje X 
graficar_freq_snp <- function(tabla_freq, titulo, subtitulo, etiquetas_x = NULL) {
  p <- tabla_freq %>%
    mutate(SNP  = factor(SNP, levels = snp_nombres),
           facet_label = paste0(SNP, " (", Gen, ")\np = ",
                                sprintf("%.3f", p_chi),
                                ifelse(p_chi < 0.05, " *", "")),
           facet_label = factor(facet_label, levels = unique(facet_label[order(SNP)]))) %>%
    pivot_longer(cols      = c(Hom_ref_pct, Het_pct, Hom_risk_pct), names_to  = "Genotipo", values_to = "Pct") %>%
    mutate(Genotipo = factor(Genotipo,
                             levels = c("Hom_ref_pct","Het_pct","Hom_risk_pct"),
                             labels = c("Hom. referencia (0)","Heterocigoto (1)","Hom. riesgo (2)"))) %>%
    ggplot(aes(x = Grupo, y = Pct, fill = Genotipo)) +
    geom_col(position = "stack", width = 0.45, alpha = 0.9) +
    geom_text(aes(label = ifelse(Pct >= 1, paste0(round(Pct), "%"), "")),
              position = position_stack(vjust = 0.5), size = 2.3, color = "black", fontface = "bold") +
    facet_wrap(~facet_label, ncol = 4) +
    scale_fill_manual(values = colores_snp, name = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%")) +
    labs(title = titulo, subtitle = subtitulo, x = NULL, y = "Porcentaje (%)") +
    tema_base +
    theme(legend.position = "bottom",
          strip.text      = element_text(size = 7, face = "bold"),
          axis.text.x     = element_text(size = 7))
  if (!is.null(etiquetas_x)) p <- p + scale_x_discrete(labels = etiquetas_x)
  p
}

# ==============================================================================
# 5.1. PREPARACIÓN DE DATOS: NORMALIZACIÓN Y VARIABLES DOSAGE
# ==============================================================================
# Aplica las funciones definidas en 5.0 al dataset para crear las variables necesarias para el análisis 
# (genotipos normalizados + dosage por SNP).

# ── 5.1.1. Normalización de genotipos en datos y datos_fh ─────────────────────
for (v in snp_raw) {
  datos_fh[[v]] <- normalizar_geno(datos_fh[[v]])
  datos[[v]]    <- normalizar_geno(datos[[v]])
}

# ── 5.1.2. Creación de variables dosage (0/1/2) ───────────────────────────────
for (i in seq_along(snp_raw)) {
  dvar             <- paste0(snp_nombres[i], "_dose")
  datos_fh[[dvar]] <- crear_dosage(snp_raw[i], snp_risk_alleles[i], datos_fh)
  datos[[dvar]]    <- crear_dosage(snp_raw[i], snp_risk_alleles[i], datos)
}
dose_vars <- paste0(snp_nombres, "_dose")

# Recrear subsets de sexo para que hereden las columnas dosage
datos_m <- filter(datos, Sexo == "Mujer")
datos_h <- filter(datos, Sexo == "Hombre")

# ── 5.1.3. Verificación ───────────────────────────────────────────────────────
cat("=== Distribución de dosage: Cohorte FH ===\n")
for (i in seq_along(snp_nombres)) {
  cat(snp_nombres[i], "(", snp_genes[i], "):",
      "N=",    sum(!is.na(datos_fh[[dose_vars[i]]])),
      "| 0:",  sum(datos_fh[[dose_vars[i]]] == 0, na.rm=TRUE),
      "| 1:",  sum(datos_fh[[dose_vars[i]]] == 1, na.rm=TRUE),
      "| 2:",  sum(datos_fh[[dose_vars[i]]] == 2, na.rm=TRUE),
      "| NA:", sum(is.na(datos_fh[[dose_vars[i]]])), "\n")
}

# ==============================================================================
# 5.2. DISTRIBUCIÓN GENOTÍPICA: 12 SNPs POR GRUPO DE ESTUDIO Y SEXO
# ==============================================================================
# Análisis descriptivo de las frecuencias genotípicas de los 12 SNPs en los tres grupos de estudio (No FH, FH Sin Evento, FH Con Evento)
#
# Preguntas:
#   - ¿Difieren las frecuencias genotípicas entre grupos? -> Test chi-cuadrado por SNP
#   - ¿Hay diferencias entre sexos en la distribución de los SNPs? -> Chi-cuadrado x sexo
#
# Nota metodológica: las 4 variantes ALOX5AP (rs10507391, rs17222842, rs9315050, rs17216473) constituyen un haplotipo 
# (haplotipo B), por lo que presentan RAF elevadas en población española y están en desequilibrio de ligamiento parcial.
# Se analizan como variantes individuales para mantener coherencia con el score original.

# ── 5.2.1. Frecuencias por grupo de estudio ──────────────────────────────────
tabla_freq_grupos <- calcular_freq_snp( datos, "Grupo_estudio", c("No FH", "FH No Evento", "FH Evento"))

cat("=== Frecuencias por grupo ===\n")
print(tabla_freq_grupos[, c("SNP","Gen","Grupo","N","Hom_ref_pct", "Het_pct","Hom_risk_pct","RAF","p_chi")],
      row.names = FALSE)

# ── 5.2.2. Tabla gt: distribución genotípica por grupo ───────────────────────
tabla_freq_grupos %>%
  mutate(Grp = case_when(
    Grupo == "No FH"        ~ "nofh",
    Grupo == "FH No Evento" ~ "sine",
    Grupo == "FH Evento"    ~ "cone")) %>%
  select(SNP, Gen, Grp, Hom_ref_pct, Het_pct, Hom_risk_pct, p_chi) %>%
  pivot_wider(
    id_cols     = c(SNP, Gen, p_chi),
    names_from  = Grp,
    values_from = c(Hom_ref_pct, Het_pct, Hom_risk_pct)) %>%
  left_join(
    tabla_freq_grupos %>%
      filter(Grupo == "FH Evento") %>%
      select(SNP, RAF),
    by = "SNP") %>%
  mutate(p_sig = p_chi < 0.05) %>%
  select(SNP, Gen, RAF,
         Hom_ref_pct_nofh, Het_pct_nofh, Hom_risk_pct_nofh,
         Hom_ref_pct_sine, Het_pct_sine, Hom_risk_pct_sine,
         Hom_ref_pct_cone, Het_pct_cone, Hom_risk_pct_cone,
         p_chi, p_sig) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 17. Distribución genotípica de los 12 SNPs por grupo de estudio**",
    subtitulo = paste0("*Cohorte completa (N=", nrow(datos), ") · No FH: n=113 · FH Sin Evento: n=609 · FH Con Evento: n=477*")) %>%
  cols_hide("p_sig") %>%
  tab_spanner(label   = "No FH (n=113)", columns = c(Hom_ref_pct_nofh, Het_pct_nofh, Hom_risk_pct_nofh)) %>%
  tab_spanner(label   = "FH Sin Evento (n=609)", columns = c(Hom_ref_pct_sine, Het_pct_sine, Hom_risk_pct_sine)) %>%
  tab_spanner(label   = "FH Con Evento (n=477)", columns = c(Hom_ref_pct_cone, Het_pct_cone, Hom_risk_pct_cone)) %>%
  cols_label(
    SNP  = "SNP",  Gen = "Gen",  RAF = "RAF",
    Hom_ref_pct_nofh  = "Hom.ref %", Het_pct_nofh  = "Het %",
    Hom_risk_pct_nofh = "Hom.R %",
    Hom_ref_pct_sine  = "Hom.ref %", Het_pct_sine  = "Het %",
    Hom_risk_pct_sine = "Hom.R %",
    Hom_ref_pct_cone  = "Hom.ref %", Het_pct_cone  = "Het %",
    Hom_risk_pct_cone = "Hom.R %",
    p_chi = md("p (χ²)")) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = p_chi, rows = p_sig)) %>%
  tab_footnote(
    footnote  = md("RAF: Risk Allele Frequency calculada en la cohorte FH Con Evento.<br>
                    χ² con simulación Monte Carlo (B=2.000).<br>
                    Hom.ref: homocigoto referencia (0 copias); Het: heterocigoto (1 copia);
                    Hom.R: homocigoto riesgo (2 copias)."),
    locations = cells_column_labels(columns = RAF)) %>%
  guardar_gt("Tabla17_Freq_SNPs_GrupoEstudio.html")

# ── 5.2.3. Figura por grupo ───────────────────────────────────────────────────
guardar_figura(
  graficar_freq_snp(
    tabla_freq_grupos,
    titulo      = "Distribución genotípica de los 12 SNPs por grupo de estudio",
    subtitulo   = "N=1.199 · * p<0,05 (χ² Monte Carlo B=2.000)",
    etiquetas_x = c("No FH"        = "No FH",
                    "FH No Evento" = "FH Sin\nEvento",
                    "FH Evento"    = "FH Con\nEvento")),
  "BarChart_SNPs_Genotipo_GrupoEstudio.png", ancho = 12, alto = 9)

# ── 5.2.4. Frecuencias por sexo (dentro de la cohorte FH) ─────────────────────
tabla_freq_sexo <- calcular_freq_snp(datos_fh, "Sexo", c("Mujer", "Hombre"))

cat("\n=== Frecuencias por sexo ===\n")
print(tabla_freq_sexo[, c("SNP","Gen","Grupo","N","Hom_ref_pct", "Het_pct","Hom_risk_pct","RAF","p_chi")],
      row.names = FALSE)

guardar_figura(
  graficar_freq_snp(
    tabla_freq_sexo,
    titulo    = "Distribución genotípica de los 12 SNPs por sexo: Cohorte FH",
    subtitulo = paste0("Mujeres n=", nrow(datos_fh_m),
                       " · Hombres n=", nrow(datos_fh_h),
                       " · * p<0,05 (χ² Monte Carlo B=2.000)")),
  "BarChart_SNPs_Genotipo_Sexo.png", ancho = 12, alto = 9)

# ── 5.2.5. Frecuencias por grupo de estudio estratificadas por sexo ───────────
# Analiza si el efecto de algún SNP entre grupos (No FH / FH Sin Evento / FH Con Evento) es específico de un sexo
# Usa datos_m y datos_h.
tabla_freq_m_grupos <- calcular_freq_snp(datos_m, "Grupo_estudio", c("No FH", "FH No Evento", "FH Evento"))
tabla_freq_h_grupos <- calcular_freq_snp(datos_h, "Grupo_estudio", c("No FH", "FH No Evento", "FH Evento"))

cat("\n=== p-valores por grupo:  Mujeres ===\n")
print(tabla_freq_m_grupos %>% select(SNP, Gen, p_chi) %>% distinct(), row.names = FALSE)

cat("\n=== p-valores por grupo:  Hombres ===\n")
print(tabla_freq_h_grupos %>% select(SNP, Gen, p_chi) %>% distinct(), row.names = FALSE)

guardar_figura(
  graficar_freq_snp(
    tabla_freq_m_grupos,
    titulo      = "Distribución genotípica de los 12 SNPs por grupo: Mujeres",
    subtitulo   = paste0(
      "No FH n=", sum(datos_m$Grupo == "No FH", na.rm=TRUE),
      " · FH Sin Evento n=", nrow(datos_fh_m),
      " · FH Con Evento n=", sum(datos_fh_m$ECV_bin == 1),
      " · * p<0,05 (χ² Monte Carlo B=2.000)"),
    etiquetas_x = c("No FH"        = "No FH",
                    "FH No Evento" = "FH Sin\nEvento",
                    "FH Evento"    = "FH Con\nEvento")),
  "BarChart_SNPs_Genotipo_Mujeres.png", ancho = 12, alto = 9)

guardar_figura(
  graficar_freq_snp(
    tabla_freq_h_grupos,
    titulo      = "Distribución genotípica de los 12 SNPs por grupo: Hombres",
    subtitulo   = paste0(
      "No FH n=", sum(datos_h$Grupo == "No FH", na.rm=TRUE),
      " · FH Sin Evento n=", nrow(datos_fh_h),
      " · FH Con Evento n=", sum(datos_fh_h$ECV_bin == 1),
      " · * p<0,05 (χ² Monte Carlo B=2.000)"),
    etiquetas_x = c("No FH"        = "No FH",
                    "FH No Evento" = "FH Sin\nEvento",
                    "FH Evento"    = "FH Con\nEvento")),
  "BarChart_SNPs_Genotipo_Hombres.png", ancho = 12, alto = 9)

# Frecuencias detalladas de los 3 SNPs significativos por grupo × sexo
snps_sig <- c("rs10455872", "rs12526453", "rs10507391")

cat("=== Frecuencias detalladas SNPs significativos × sexo ===\n")
bind_rows(
  tabla_freq_m_grupos %>% mutate(Sexo = "Mujeres"),
  tabla_freq_h_grupos %>% mutate(Sexo = "Hombres")) %>%
  filter(SNP %in% snps_sig) %>%
  select(Sexo, SNP, Gen, Grupo, N, Hom_ref_pct, Het_pct, Hom_risk_pct, RAF, p_chi) %>%
  arrange(SNP, Sexo, Grupo) %>%
  print(row.names = FALSE)

# Figura resumen SNPs significativos
# Panel A: rs10455872 (LPA) - Global
p_lpa_global <- graficar_freq_snp(
  tabla_freq_grupos %>% filter(SNP == "rs10455872"),
  titulo = "A. Población Global", subtitulo = NULL,
  etiquetas_x = c("No FH" = "No FH", "FH No Evento" = "FH Sin\nEvento", "FH Evento" = "FH Con\nEvento"))

# Panel B: rs10455872 (LPA) - Mujeres
p_lpa_muj <- graficar_freq_snp(
  tabla_freq_m_grupos %>% filter(SNP == "rs10455872"),
  titulo = "B. Mujeres", subtitulo = NULL,
  etiquetas_x = c("No FH" = "No FH", "FH No Evento" = "FH Sin\nEvento", "FH Evento" = "FH Con\nEvento"))

# Panel C: rs12526453 (PHACTR1) - Mujeres
p_phactr1_muj <- graficar_freq_snp(
  tabla_freq_m_grupos %>% filter(SNP == "rs12526453"),
  titulo = "C. Mujeres", subtitulo = NULL,
  etiquetas_x = c("No FH" = "No FH", "FH No Evento" = "FH Sin\nEvento", "FH Evento" = "FH Con\nEvento"))

# Panel D: rs10507391 (ALOX5AP) - Hombres
p_alox_hom <- graficar_freq_snp(
  tabla_freq_h_grupos %>% filter(SNP == "rs10507391"),
  titulo = "D. Hombres", subtitulo = NULL,
  etiquetas_x = c("No FH" = "No FH", "FH No Evento" = "FH Sin\nEvento", "FH Evento" = "FH Con\nEvento"))

# Uniir los 4 gráficos y unificar la leyenda
figura_resumen <- (p_lpa_global | p_lpa_muj) / (p_phactr1_muj | p_alox_hom) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
guardar_figura(figura_resumen, "Panel_SNPs_Significativos.png", ancho = 10, alto = 8)


# ==============================================================================
# 5.3. CARGA ALÉLICA TOTAL
# ==============================================================================
# Suma de los 12 dosages individuales (rango 0 - 24), sin ponderar por los pesos poblacionales del Score. Mide el 
# conteo bruto de alelos de riesgo por individuo, tratando todos los SNPs como igualmente informativos.
#
# Hipótesis: si el burden total difiere entre FH Sin Evento y FH Con Evento, los SNPs tienen valor discriminante en
# FH con independencia de los pesos del score. 

# ── 5.3.1. Creación de la variable de carga alélica ───────────────────────────
datos_fh$carga_total <- rowSums(datos_fh[, dose_vars], na.rm = FALSE)

# Recrear subsets de sexo para heredar carga_total
datos_fh_m <- filter(datos_fh, Sexo == "Mujer")
datos_fh_h <- filter(datos_fh, Sexo == "Hombre")

cat("=== Distribución de la carga alélica total (FH) ===\n")
cat("Rango:", range(datos_fh$carga_total, na.rm=TRUE), "\n")
cat("Mediana (Sin Evento):", median(datos_fh$carga_total[datos_fh$Grupo_estudio=="FH No Evento"], na.rm=TRUE), "\n")
cat("Mediana (Con Evento):", median(datos_fh$carga_total[datos_fh$Grupo_estudio=="FH Evento"],    na.rm=TRUE), "\n")

# ── 5.3.2. Tests estadísticos ──────────────────────────────────────────────────
# Global FH
wilcox_carga <- wilcox.test(carga_total ~ Grupo_estudio, data = datos_fh)
pw_carga     <- pairwise.wilcox.test(datos_fh$carga_total, datos_fh$Grupo_estudio, p.adjust.method = "bonferroni")

# Por sexo
wilcox_carga_m <- wilcox.test(carga_total ~ Grupo_estudio, data = datos_fh_m)
wilcox_carga_h <- wilcox.test(carga_total ~ Grupo_estudio, data = datos_fh_h)
pw_carga_m     <- pairwise.wilcox.test(datos_fh_m$carga_total, datos_fh_m$Grupo_estudio, p.adjust.method = "bonferroni")
pw_carga_h     <- pairwise.wilcox.test(datos_fh_h$carga_total, datos_fh_h$Grupo_estudio, p.adjust.method = "bonferroni")

cat("\n=== Wilcoxon carga alélica (FH global) ===\n")
cat("p =", round(wilcox_carga$p.value, 4), "\n")
cat("Post-hoc (Wilcoxon-Bonferroni):\n")
print(wilcox_carga$p.value)

cat("\n=== Wilcoxon carga alélica: Mujeres ===\n")
cat("p =", round(wilcox_carga_m$p.value, 4), "\n")
print(wilcox_carga_m$p.value)

cat("\n=== Wilcoxon carga alélica: Hombres ===\n")
cat("p =", round(wilcox_carga_h$p.value, 4), "\n")
print(wilcox_carga_h$p.value)

cat("=== Medianas de carga alélica por grupo y sexo ===\n")
datos_fh %>%
  group_by(Sexo, Grupo_estudio) %>%
  summarise(
    N       = n(),
    Mediana = median(carga_total, na.rm=TRUE),
    Media   = round(mean(carga_total, na.rm=TRUE), 2),
    Q1      = quantile(carga_total, 0.25, na.rm=TRUE),
    Q3      = quantile(carga_total, 0.75, na.rm=TRUE),
    .groups = "drop") %>%
  print(n=Inf)


# ── 5.3.3. Figura: boxplot carga alélica × grupo × sexo ───────────────────────
y_max_carga <- max(datos_fh$carga_total, na.rm=TRUE) * 1.18

signif_carga <- data.frame(
  Sexo       = factor(c("Mujer","Hombre"), levels = c("Mujer","Hombre")),
  xmin       = c(1, 1),
  xmax       = c(2, 2),
  y_position = c(y_max_carga * 0.93, y_max_carga * 0.93),
  label      = c(sig_label(wilcox_carga_m$p.value), sig_label(wilcox_carga_h$p.value)))

fig_carga <- datos_fh %>%
  filter(!is.na(carga_total)) %>%
  ggplot(aes(x = Grupo_estudio, y = carga_total, fill = Grupo_estudio)) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.4, color = "grey30") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_signif(data         = signif_carga,
              aes(annotations = label),
              xmin         = signif_carga$xmin,
              xmax         = signif_carga$xmax,
              y_position   = signif_carga$y_position,
              manual       = TRUE, inherit.aes = FALSE,
              textsize     = 3, vjust = 0.3, size = 0.3) +
  facet_wrap(~Sexo, labeller = labeller(Sexo = c("Mujer"  = "Mujeres", "Hombre" = "Hombres"))) +
  scale_fill_manual(values = colores_ecv) +
  scale_x_discrete(labels = c("FH No Evento" = "FH Sin\nEvento", "FH Evento" = "FH Con\nEvento")) +
  scale_y_continuous(limits = c(NA, y_max_carga)) +
  labs(title    = "Carga alélica total (∑ 12 SNPs) por grupo de evento y sexo",
       subtitle = paste0("Mujeres: Wilcoxon p=", round(wilcox_carga_m$p.value, 3), "  |  Hombres: Wilcoxon p=", round(wilcox_carga_h$p.value, 3)),
       caption  = paste0("Carga alélica = suma de dosages (0/1/2) de los 12 SNPs, rango 0-24. ",
                         "Wilcoxon global FH: p=", round(wilcox_carga$p.value, 3)),
       x = NULL, y = "Carga alélica total (0-24)") +
  tema_base +
  theme(legend.position = "none", 
        strip.text = element_text(face = "plain", size = 9),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.major.y = element_line(linewidth = 0.3))

guardar_figura(fig_carga, "Boxplot_CargaAlelica_GrupoSexo.png", ancho = 5.5, alto = 4.5)


# ==============================================================================
# 5.4. GLMM POR SNP INDIVIDUAL + CORRECCIÓN FDR (BENJAMINI-HOCHBERG)
# ==============================================================================
# Para cada uno de los 12 SNPs: GLMM con las mismas covariables del modelo base (Bloque 3) + dosage del SNP como
# predictor aditivo (0/1/2 alelos de riesgo).
# 
# FDR Benjamini-Hochberg: controla la tasa de falsos positivos esperada (más adecuado que Bonferroni para análisis
# exploratorio con 12 tests correlacionados).
# 
# Análisis: global FH + estratificado por sexo (sin Sexo en las covariables).

# ── 5.4.1. Función auxiliar ────────────────────────────────────────────────────
# Ajusta 12 GLMMs (uno por SNP) sobre cualquier dataset y devuelve tabla de ORs
# sin_sexo: TRUE para análisis estratificados (elimina Sexo de covariables)
glmm_por_snp <- function(df, sin_sexo = FALSE) {
  cov <- if (sin_sexo) {
    "scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  } else {
    "Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin"
  }
  
  do.call(rbind, lapply(seq_along(snp_nombres), function(i) {
    snp_var <- dose_vars[i]
    f <- as.formula(paste("ECV_bin ~", cov, "+", snp_var, "+ (1 | ID_cluster)"))
    
    modelo <- tryCatch(
      glmer(f, data = df, family = binomial, control = ctrl),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        glmer(f, data = df, family = binomial, control = ctrl)))
    
    if (is.null(modelo)) {
      return(data.frame(SNP=snp_nombres[i], Gen=snp_genes[i], OR=NA, IC_2.5=NA, IC_97.5=NA, p_raw=NA))}
    
    coefs <- summary(modelo)$coefficients
    ci    <- tryCatch(
      exp(confint(modelo, method = "Wald", parm = "beta_")),
      error = function(e) matrix(NA, nrow=length(fixef(modelo)), ncol=2, dimnames=list(names(fixef(modelo)), c("2.5 %","97.5 %"))))
    
    data.frame(
      SNP     = snp_nombres[i],
      Gen     = snp_genes[i],
      OR      = round(exp(coefs[snp_var, "Estimate"]), 3),
      IC_2.5  = round(ci[snp_var, "2.5 %"],  3),
      IC_97.5 = round(ci[snp_var, "97.5 %"], 3),
      p_raw   = coefs[snp_var, "Pr(>|z|)"])
  }))
}

# ── 5.4.2. Ajuste de modelos ───────────────────────────────────────────────────
res_snp_global <- glmm_por_snp(datos_fh, sin_sexo = FALSE)
res_snp_m      <- glmm_por_snp(datos_fh_m, sin_sexo = TRUE)
res_snp_h      <- glmm_por_snp(datos_fh_h, sin_sexo = TRUE)

# ── 5.4.3. Corrección FDR (Benjamini-Hochberg) ────────────────────────────────
res_snp_global$p_FDR <- p.adjust(res_snp_global$p_raw, method = "BH")
res_snp_m$p_FDR      <- p.adjust(res_snp_m$p_raw,      method = "BH")
res_snp_h$p_FDR      <- p.adjust(res_snp_h$p_raw,      method = "BH")

# ── 5.4.4. Resultados ─────────────────────────────────────────────────────────
cat("\n=== GLMM por SNP — FH Global ===\n")
print(res_snp_global %>%
        mutate(p_raw = round(p_raw, 4), p_FDR = round(p_FDR, 4), sig_FDR = ifelse(p_FDR < 0.05, "*", "")) %>%
        select(SNP, Gen, OR, IC_2.5, IC_97.5, p_raw, p_FDR, sig_FDR),
      row.names = FALSE)

cat("\n=== GLMM por SNP — Mujeres ===\n")
print(res_snp_m %>%
        mutate(p_raw = round(p_raw, 4), p_FDR = round(p_FDR, 4), sig_FDR = ifelse(p_FDR < 0.05, "*", "")) %>%
        select(SNP, Gen, OR, IC_2.5, IC_97.5, p_raw, p_FDR, sig_FDR),
      row.names = FALSE)

cat("\n=== GLMM por SNP — Hombres ===\n")
print(res_snp_h %>%
        mutate(p_raw = round(p_raw, 4), p_FDR = round(p_FDR, 4), sig_FDR = ifelse(p_FDR < 0.05, "*", "")) %>%
        select(SNP, Gen, OR, IC_2.5, IC_97.5, p_raw, p_FDR, sig_FDR),
      row.names = FALSE)

# ── 5.4.5. Tabla 18: GLMM por SNP x FDR ───────────────────────────────────────
tabla18_datos <- bind_rows(
  res_snp_global %>% mutate(Analisis = "Global FH"),
  res_snp_m      %>% mutate(Analisis = "Mujeres"),
  res_snp_h      %>% mutate(Analisis = "Hombres")) %>%
  mutate(
    Analisis  = factor(Analisis, levels = c("Global FH", "Mujeres", "Hombres")),
    SNP       = factor(SNP, levels = snp_nombres),
    OR_f      = sprintf("%.3f", OR),
    IC_95     = paste0("(", sprintf("%.3f", IC_2.5), " – ", sprintf("%.3f", IC_97.5), ")"),
    p_raw_f   = ifelse(p_raw < 0.001, "<0.001", sprintf("%.3f", p_raw)),
    p_FDR_f   = ifelse(p_FDR < 0.001, "<0.001", sprintf("%.3f", p_FDR)),
    p_raw_sig = p_raw < 0.05,
    p_FDR_sig = p_FDR < 0.05) %>%
  arrange(Analisis, SNP)

cat("=== Tabla 18 ===\n")
print(tabla18_datos %>% select(Analisis, SNP, Gen, OR, IC_95, p_raw_f, p_FDR_f))

tabla18_datos %>%
  select(Analisis, SNP, Gen, OR_f, IC_95, p_raw_f, p_FDR_f, p_raw_sig, p_FDR_sig) %>%
  gt(groupname_col = "Analisis") %>%
  gt_estilo(
    titulo    = "**Tabla 18. GLMM por SNP individual con corrección FDR**",
    subtitulo = paste0(
      "*OR por cada alelo de riesgo adicional (dosage 0/1/2), ajustado por el modelo base*<br>",
      "*Global FH N=", nrow(datos_fh),
      " · Mujeres N=", nrow(datos_fh_m),
      " · Hombres N=", nrow(datos_fh_h), "*")) %>%
  cols_label(
    SNP     = "SNP",
    Gen     = "Gen",
    OR_f    = "OR",
    IC_95   = "IC 95%",
    p_raw_f = "p (crudo)",
    p_FDR_f = md("p (FDR-BH)")) %>%
  cols_align(align = "left",  columns = c(SNP, Gen, IC_95)) %>%
  cols_align(align = "right", columns = c(OR_f, p_raw_f, p_FDR_f)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_raw_f", rows = p_raw_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = "p_FDR_f", rows = p_FDR_sig)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote  = md("OR ajustados por Edad, cLDL, cHDL, Lp(a), HTA y DM (IC 95% Wald).<br>
                    Mujeres y Hombres: Sexo excluido como covariable (análisis estratificado).<br>
                    FDR: corrección Benjamini-Hochberg para 12 comparaciones por análisis. 
                    Efecto aleatorio: (1 | ID_cluster)."),
    locations = cells_column_labels(columns = "OR_f")) %>%
  cols_hide(c("p_raw_sig", "p_FDR_sig")) %>%
  guardar_gt("Tabla18_GLMM_SNP_FDR.html")


# ==============================================================================
# 5.5. RECALIBRACIÓN FH - ESPECÍFICA DEL SCORE
# ==============================================================================
# En lugar de usar los pesos del score (calibrados en población general), se ajusta un modelo con los 
# 12 SNPs simultáneamente, dejando que la propia cohorte FH estime los pesos relativos. Se valida por cross-validation
# para evitar sobreajuste.

# ── 5.5.1. Comprobación de colinealidad (VIF) ─────────────────────────────────
# Las 4 variantes ALOX5AP forman un haplotipo y pueden estar correlacionadas.
# VIF > 5 indica colinealidad problemática; VIF > 10 es severa.
modelo_vif <- lm(as.formula(paste("ECV_bin ~", paste(dose_vars, collapse = " + "))), data = datos_fh)

cat("=== VIF: Colinealidad entre los 12 SNPs ===\n")
print(round(vif(modelo_vif), 2))

# ── 5.5.2. Modelo conjunto: los 12 SNPs simultáneos ───────────────────────────
# Mismas covariables del modelo base (Bloque 3). Todos los SNPs entran juntos como predictores aditivos: sus 
# coeficientes son los pesos FH-específicos.

f_conjunto_global <- as.formula(paste(
  "ECV_bin ~ Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +",
  "HTA_bin + DM_bin +", paste(dose_vars, collapse = " + "), "+ (1 | ID_cluster)"))

f_conjunto_sexo <- as.formula(paste(
  "ECV_bin ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) +",
  "HTA_bin + DM_bin +", paste(dose_vars, collapse = " + "), "+ (1 | ID_cluster)"))

modelo_conjunto_global <- glmer(f_conjunto_global, data = datos_fh, family = binomial, control = ctrl)
modelo_conjunto_m      <- glmer(f_conjunto_sexo, data = datos_fh_m, family = binomial, control = ctrl)
modelo_conjunto_h      <- glmer(f_conjunto_sexo, data = datos_fh_h, family = binomial, control = ctrl)

cat("=== Modelo conjunto: FH Global ===\n")
print(tabla_ors(modelo_conjunto_global))
cat("\n=== Modelo conjunto: Mujeres ===\n")
print(tabla_ors(modelo_conjunto_m))
cat("\n=== Modelo conjunto: Hombres ===\n")
print(tabla_ors(modelo_conjunto_h))

# ── 5.5.4. Validación cruzada k-fold respetando estructura familiar ───────────
# ── 5.5.4.1. Asignación de folds ──────────────────────────────────────────────
# Los folds se asignan a nivel de ID_cluster (familia), no de paciente individual, para evitar fuga de información 
# cuando varios miembros de una familia caen en folds distintos. 
asignar_folds_familia <- function(df, k = 10) {
  familias <- unique(df$ID_cluster)
  fold_familia <- sample(rep(1:k, length.out = length(familias)))
  names(fold_familia) <- as.character(familias)
  fold_familia[as.character(df$ID_cluster)]
}

# Verificación: distribución de pacientes por fold (global FH)
folds_global <- asignar_folds_familia(datos_fh, k = 10)
cat("=== Distribución de pacientes por fold: FH Global ===\n")
print(table(folds_global))

cat("\n=== Verificación: todas las familias en un único fold ===\n")
chequeo <- datos_fh %>%
  mutate(fold = folds_global) %>%
  group_by(ID_cluster) %>%
  summarise(n_folds_distintos = n_distinct(fold)) %>%
  pull(n_folds_distintos)
cat("Máximo de folds distintos por familia (debe ser 1):", max(chequeo), "\n")

# ── 5.5.4.2. Función de validación cruzada k-fold con efecto aleatorio ───────
# Ajusta el modelo en k-1 folds y predice sobre el fold reservado, repitiendo k veces. Devuelve el vector de 
# predicciones out-of-fold y el AUC resultante. Los modelos que no convergen en un fold se omiten con aviso (tryCatch).
cv_glmm_auc <- function(formula, df, folds, k = 10) {
  pred_oof <- rep(NA_real_, nrow(df))
  n_fallos <- 0
  
  for (i in 1:k) {
    train <- df[folds != i, ]
    test  <- df[folds == i, ]
    
    modelo <- tryCatch(
      suppressWarnings(glmer(formula, data = train, family = binomial, control = ctrl)),
      error = function(e) NULL)
    
    if (is.null(modelo)) {
      n_fallos <- n_fallos + 1
      next
    }
    
    # predict con allow.new.levels=TRUE: ID_cluster del fold de test no estaba en el entrenamiento
    pred_oof[folds == i] <- predict(modelo, newdata = test, type = "response", allow.new.levels = TRUE)
  }
  
  if (n_fallos > 0) cat("  (", n_fallos, "de", k, "folds no convergieron)\n")
  
  auc_val <- as.numeric(auc(roc(df$ECV_bin[!is.na(pred_oof)], pred_oof[!is.na(pred_oof)], quiet = TRUE)))
  list(pred = pred_oof, auc = auc_val, n_fallos = n_fallos)
}

# ── 5.5.4.3. Folds para mujeres y hombres ─────────────────────────────────────
folds_m <- asignar_folds_familia(datos_fh_m, k = 10)
folds_h <- asignar_folds_familia(datos_fh_h, k = 10)

# ── 5.5.4.4. Fórmulas de los 3 modelos a comparar ─────────────────────────────
# f_conjunto_global y f_conjunto_sexo ya definidos en 5.5.3 
f_base_global <- as.formula(
  "ECV_bin ~ Sexo + scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin + (1 | ID_cluster)")
f_grs_global  <- update(f_base_global, . ~ . + scale(GRS))

f_base_sexo <- as.formula(
  "ECV_bin ~ scale(Edad_inclusion) + scale(cLDL_0) + scale(cHDL_0) + scale(LpA_0) + HTA_bin + DM_bin + (1 | ID_cluster)")
f_grs_sexo  <- update(f_base_sexo, . ~ . + scale(GRS))

# ── 5.5.4.5. Ejecución de la validación cruzada ───────────────────────────────
cat("=== CV: FH Global ===\n")
cv_base_global     <- cv_glmm_auc(f_base_global,     datos_fh, folds_global)
cv_grs_global      <- cv_glmm_auc(f_grs_global,      datos_fh, folds_global)
cv_conjunto_global <- cv_glmm_auc(f_conjunto_global, datos_fh, folds_global)

cat("\n=== CV: Mujeres ===\n")
cv_base_m     <- cv_glmm_auc(f_base_sexo,     datos_fh_m, folds_m)
cv_grs_m      <- cv_glmm_auc(f_grs_sexo,      datos_fh_m, folds_m)
cv_conjunto_m <- cv_glmm_auc(f_conjunto_sexo, datos_fh_m, folds_m)

cat("\n=== CV: Hombres ===\n")
cv_base_h     <- cv_glmm_auc(f_base_sexo,     datos_fh_h, folds_h)
cv_grs_h      <- cv_glmm_auc(f_grs_sexo,      datos_fh_h, folds_h)
cv_conjunto_h <- cv_glmm_auc(f_conjunto_sexo, datos_fh_h, folds_h)

cat("\n=== AUC out-of-fold (resumen) ===\n")
cat("Global: Base:", round(cv_base_global$auc,4),
    "| GRS:", round(cv_grs_global$auc,4),
    "| Conjunto:", round(cv_conjunto_global$auc,4), "\n")
cat("Mujeres: Base:", round(cv_base_m$auc,4),
    "| GRS:", round(cv_grs_m$auc,4),
    "| Conjunto:", round(cv_conjunto_m$auc,4), "\n")
cat("Hombres: Base:", round(cv_base_h$auc,4),
    "| GRS:", round(cv_grs_h$auc,4),
    "| Conjunto:", round(cv_conjunto_h$auc,4), "\n")

# ── 5.5.4.6. Repetición de la CV en hombres (estabilidad de la estimación) ────
# Repite el 10-fold CV con 5 particiones distintas (semillas distintas) para comprobar si la mejora de AUC 
# observada es estable o depende de la partición concreta.
repetir_cv_hombres <- function(n_rep = 5) {
  resultados <- data.frame()
  for (r in 1:n_rep) {
    set.seed(2026 + r)
    folds_h_r <- asignar_folds_familia(datos_fh_h, k = 10)
    
    auc_base_r     <- cv_glmm_auc(f_base_sexo,     datos_fh_h, folds_h_r)$auc
    auc_grs_r      <- cv_glmm_auc(f_grs_sexo,      datos_fh_h, folds_h_r)$auc
    auc_conjunto_r <- cv_glmm_auc(f_conjunto_sexo, datos_fh_h, folds_h_r)$auc
    
    resultados <- rbind(resultados, data.frame(
      Repeticion = r, Base = auc_base_r, GRS = auc_grs_r, Conjunto = auc_conjunto_r))
  }
  resultados
}

cv_rep_hombres <- repetir_cv_hombres(n_rep = 5)

cat("=== CV repetida (5 × 10-fold) Hombres ===\n")
print(cv_rep_hombres)

cat("\n=== Resumen (media ± DE) ===\n")
cv_rep_hombres %>% summarise(across(c(Base, GRS, Conjunto), list(media = mean, sd = sd))) %>% print()

# ── 5.5.5. Tabla 19: AUC validación cruzada: comparación de modelos ──────────
tabla19_datos <- data.frame(
  Modelo   = c("Base (sin SNPs)", "GRS original (pesos poblacionales)", "Recalibrado (12 SNPs FH-específicos)"),
  AUC_G  = c(cv_base_global$auc, cv_grs_global$auc, cv_conjunto_global$auc),    # Global FH
  AUC_M  = c(cv_base_m$auc, cv_grs_m$auc, cv_conjunto_m$auc),                   # Mujeres
  AUC_H  = c(mean(cv_rep_hombres$Base), mean(cv_rep_hombres$GRS),               # Hombres: media de 5 repeticiones
             mean(cv_rep_hombres$Conjunto)),
  AUC_H_sd = c(sd(cv_rep_hombres$Base), sd(cv_rep_hombres$GRS), sd(cv_rep_hombres$Conjunto))) %>%
  mutate(
    dAUC_G  = round(AUC_G - AUC_G[1], 4),
    dAUC_M  = round(AUC_M - AUC_M[1], 4),
    dAUC_H  = round(AUC_H - AUC_H[1], 4),
    AUC_G_f = sprintf("%.3f", AUC_G),
    AUC_M_f = sprintf("%.3f", AUC_M),
    AUC_H_f = paste0(sprintf("%.3f", AUC_H), " ± ", sprintf("%.3f", AUC_H_sd)),
    dAUC_G_f = ifelse(dAUC_G == 0, "—", ifelse(dAUC_G > 0, paste0("+", sprintf("%.3f", dAUC_G)), sprintf("%.3f", dAUC_G))),
    dAUC_M_f = ifelse(dAUC_M == 0, "—", ifelse(dAUC_M > 0, paste0("+", sprintf("%.3f", dAUC_M)), sprintf("%.3f", dAUC_M))),
    dAUC_H_f = ifelse(dAUC_H == 0, "—", ifelse(dAUC_H > 0, paste0("+", sprintf("%.3f", dAUC_H)), sprintf("%.3f", dAUC_H))),
    fila_base = Modelo == "Base (sin SNPs)")

tabla19_datos %>%
  select(Modelo, AUC_G_f, dAUC_G_f, AUC_M_f, dAUC_M_f, AUC_H_f, dAUC_H_f, fila_base) %>%
  gt() %>%
  gt_estilo(
    titulo    = "**Tabla 19. Capacidad predictiva: validación cruzada 10-fold**",
    subtitulo = paste0(
      "*Comparación del AUC out-of-fold entre modelo base, GRS original y score recalibrado*<br>",
      "*Global FH N=", nrow(datos_fh), " · Mujeres N=", nrow(datos_fh_m), " · Hombres N=", nrow(datos_fh_h), "*")) %>%
  tab_spanner(label   = paste0("Global FH (N=", nrow(datos_fh), ")"), columns = c(AUC_G_f, dAUC_G_f)) %>%
  tab_spanner(label   = paste0("Mujeres (N=", nrow(datos_fh_m), ")"), columns = c(AUC_M_f, dAUC_M_f)) %>%
  tab_spanner(label   = paste0("Hombres (N=", nrow(datos_fh_h), ")"), columns = c(AUC_H_f, dAUC_H_f)) %>%
  cols_label(
    Modelo   = "Modelo",
    AUC_G_f  = "AUC",  dAUC_G_f = "ΔAUC",
    AUC_M_f  = "AUC",  dAUC_M_f = "ΔAUC",
    AUC_H_f  = md("AUC (media ± DE)"), dAUC_H_f = "ΔAUC") %>%
  cols_hide("fila_base") %>%
  cols_align(align = "left",  columns = Modelo) %>%
  cols_align(align = "right", columns = c(AUC_G_f, dAUC_G_f, AUC_M_f, dAUC_M_f, AUC_H_f, dAUC_H_f)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = fila_base)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_footnote(
    footnote = md("AUC estimado mediante validación cruzada 10-fold con folds asignados a nivel de familia (ID_cluster)
                   para respetar la estructura de dependencia intrafamiliar.<br>
                   ΔAUC = diferencia respecto al modelo base. Hombres: media ± DE de 5 repeticiones con particiones distintas.<br>
                   Modelo base: Edad, cLDL, cHDL, Lp(a), HTA, DM (± Sexo en análisis global).<br>
                   Recalibrado: 12 SNPs del CARDIO inCode-Score con pesos re-estimados en la cohorte FH."),
    locations = cells_column_labels(columns = "AUC_G_f")) %>%
  guardar_gt("Tabla19_CV_AUC_Recalibracion.html")



