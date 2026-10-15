packages <- c(
  "openxlsx",
  "dplyr",
  "tidyr",
  "ggplot2"
)

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(char = packages)

# ============================================================
# CONFIGURAÇÕES
# ============================================================

arquivo_excel <- file.path(
  "Resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes",
  "Todos_os_modelos",
  "resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes.xlsx"
)

aba_dados <- "Scenario_Info"

if (!file.exists(arquivo_excel)) {
  stop(
    paste0(
      "Arquivo Excel não encontrado em:\n",
      arquivo_excel,
      "\n\nConfira o caminho em 'arquivo_excel'."
    )
  )
}

abas_excel <- openxlsx::getSheetNames(arquivo_excel)

if (!(aba_dados %in% abas_excel)) {
  stop(
    paste0(
      "A aba '", aba_dados, "' não foi encontrada no arquivo Excel.\n",
      "Abas disponíveis:\n",
      paste(abas_excel, collapse = ", ")
    )
  )
}

# ============================================================
# FUNÇÃO AUXILIAR
# ============================================================

to_numeric_safe <- function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
}

# ============================================================
# LENDO OS DADOS
# ============================================================

scenario_info <- openxlsx::read.xlsx(
  xlsxFile = arquivo_excel,
  sheet = aba_dados,
  check.names = FALSE
)

# ============================================================
# CHECANDO COLUNAS NECESSÁRIAS
# ============================================================

required_cols <- c(
  "TargetInbreeding",
  "ScenarioCandidateN",
  "ScenarioCandidateMeanF",
  "ScenarioCandidateSDF",
  "ScenarioCandidateMeanPheno",
  "ScenarioCandidateSDPheno"
)

missing_cols <- setdiff(required_cols, names(scenario_info))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "As seguintes colunas não foram encontradas na aba 'Scenario_Info': ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# PREPARANDO RESUMO DOS CANDIDATOS
# ============================================================

resumo_candidatos <- scenario_info %>%
  mutate(
    TargetInbreeding = round(to_numeric_safe(TargetInbreeding), 2),
    ScenarioCandidateN = as.integer(ScenarioCandidateN),
    ScenarioLabel = paste0("F = ", sprintf("%.2f", TargetInbreeding)),
    
    MeanF = to_numeric_safe(ScenarioCandidateMeanF),
    SDF = to_numeric_safe(ScenarioCandidateSDF),
    
    MeanPheno = to_numeric_safe(ScenarioCandidateMeanPheno),
    SDPheno = to_numeric_safe(ScenarioCandidateSDPheno)
  ) %>%
  select(
    TargetInbreeding,
    ScenarioLabel,
    ScenarioCandidateN,
    MeanF,
    SDF,
    MeanPheno,
    SDPheno
  ) %>%
  arrange(TargetInbreeding)

print(dplyr::as_tibble(resumo_candidatos), n = Inf)

# ============================================================
# ORGANIZANDO PARA O GRÁFICO
# ============================================================

dados_plot <- resumo_candidatos %>%
  pivot_longer(
    cols = c(
      MeanF,
      SDF,
      MeanPheno,
      SDPheno
    ),
    names_to = c(".value", "Trait"),
    names_pattern = "^(Mean|SD)(.*)$"
  ) %>%
  mutate(
    Trait = dplyr::case_when(
      Trait == "F" ~ "Genomic endogamy",
      Trait == "Pheno" ~ "Phenotypic value",
      TRUE ~ Trait
    ),
    Trait = factor(
      Trait,
      levels = c(
        "Genomic endogamy",
        "Phenotypic value"
      )
    ),
    ymin = Mean - SD,
    ymax = Mean + SD
  )

# ============================================================
# GRÁFICO: MÉDIA ± DESVIO-PADRÃO
# ============================================================

p_traits_candidatos <- ggplot(
  dados_plot,
  aes(
    x = ScenarioLabel,
    y = Mean,
    fill = ScenarioLabel
  )
) +
  geom_col(
    width = 0.70,
    color = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax
    ),
    width = 0.22,
    linewidth = 0.75
  ) +
  facet_wrap(
    ~ Trait,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    x = "Target inbreeding class",
    y = "Mean ± SD",
    fill = "Target inbreeding",
    title = "Mean traits of scenario candidate individuals by inbreeding class",
    subtitle = "Bars represent means of the scenario candidate individuals; error bars represent ±1 SD"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "none"
  )

print(p_traits_candidatos)

# ============================================================
# SALVANDO GRÁFICO E RESUMO
# ============================================================

pasta_saida <- file.path(
  dirname(arquivo_excel),
  "Graficos_Candidatos_Cenario"
)

if (!dir.exists(pasta_saida)) {
  dir.create(pasta_saida, recursive = TRUE)
}

arquivo_grafico <- file.path(
  pasta_saida,
  "media_endogamia_e_fenotipo_candidatos_por_faixa_endogamia.png"
)

ggsave(
  filename = arquivo_grafico,
  plot = p_traits_candidatos,
  width = 10,
  height = 5,
  dpi = 500,
  bg = "white"
)

arquivo_resumo <- file.path(
  pasta_saida,
  "resumo_endogamia_e_fenotipo_candidatos_por_faixa_endogamia.xlsx"
)

openxlsx::write.xlsx(
  list(
    Resumo_Candidatos = resumo_candidatos,
    Dados_Plot = dados_plot
  ),
  file = arquivo_resumo,
  overwrite = TRUE
)

cat("\nGráfico salvo em:\n")
cat(normalizePath(arquivo_grafico), "\n")

cat("\nResumo salvo em:\n")
cat(normalizePath(arquivo_resumo), "\n")
