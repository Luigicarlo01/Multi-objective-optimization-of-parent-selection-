packages <- c(
  "openxlsx",
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr"
)

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(char = packages)

# ============================================================
# CONFIGURAÇÕES
# ============================================================

# Use "all" para considerar todas as gerações/ciclos
# Use um número, por exemplo 5, para usar apenas uma geração
generation_to_use <- "all"
# generation_to_use <- 5

pasta_saida_base <- "Resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes"

arquivo_excel <- file.path(
  pasta_saida_base,
  "Todos_os_modelos",
  "resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes.xlsx"
)

if (!file.exists(arquivo_excel)) {
  stop(
    paste0(
      "Não encontrei o arquivo Excel em:\n",
      arquivo_excel,
      "\n\nConfira se o caminho está correto."
    )
  )
}

abas_excel <- openxlsx::getSheetNames(arquivo_excel)

if (!("UpSet_Input" %in% abas_excel)) {
  stop(
    paste0(
      "Não encontrei a aba 'UpSet_Input' no arquivo Excel.\n",
      "Abas disponíveis:\n",
      paste(abas_excel, collapse = ", ")
    )
  )
}

# ============================================================
# LENDO A BASE
# ============================================================

upset_input_xlsx <- openxlsx::read.xlsx(
  xlsxFile = arquivo_excel,
  sheet = "UpSet_Input"
)

# ============================================================
# FUNÇÃO PARA NOMES BONITOS
# ============================================================

pretty_model_name <- function(x) {
  x <- as.character(x)
  
  dplyr::case_when(
    grepl("^Strategy[_ ]1", x) ~ "Strategy 1",
    grepl("^Strategy[_ ]2", x) ~ "Strategy 2",
    grepl("^SubsetOCS", x) ~ "SubsetOCS",
    grepl("^EA", x) ~ "EA",
    grepl("^Local", x) ~ "Local",
    grepl("^SA", x) ~ "SA",
    grepl("^GA", x) ~ "GA",
    grepl("^GRASP", x) ~ "GRASP",
    grepl("^TABU", x) ~ "TABU",
    grepl("^EDA", x) ~ "EDA",
    grepl("^PSO", x) ~ "PSO",
    TRUE ~ x
  )
}

# ============================================================
# CHECAGENS E PREPARO
# ============================================================

if (!("Replication" %in% names(upset_input_xlsx))) {
  upset_input_xlsx <- upset_input_xlsx %>%
    mutate(Replication = 1L)
}

metadata_cols <- c(
  "Replication",
  "TargetInbreeding",
  "ObservedInitialInbreeding",
  "ScenarioCycleFound",
  "Generation",
  "IndividualID"
)

required_cols <- c(
  "Replication",
  "TargetInbreeding",
  "IndividualID"
)

missing_required <- setdiff(required_cols, names(upset_input_xlsx))

if (length(missing_required) > 0) {
  stop(
    paste0(
      "As seguintes colunas obrigatórias não foram encontradas na aba 'UpSet_Input': ",
      paste(missing_required, collapse = ", ")
    )
  )
}

set_cols <- setdiff(names(upset_input_xlsx), metadata_cols)

if (length(set_cols) < 2) {
  stop("Menos de dois modelos foram encontrados na aba 'UpSet_Input'.")
}

model_order <- set_cols

upset_input_xlsx <- upset_input_xlsx %>%
  mutate(
    TargetInbreeding = round(as.numeric(TargetInbreeding), 2),
    Replication = as.integer(Replication)
  )

if ("Generation" %in% names(upset_input_xlsx)) {
  upset_input_xlsx <- upset_input_xlsx %>%
    mutate(Generation = as.integer(Generation))
}

upset_input_xlsx <- upset_input_xlsx %>%
  mutate(
    across(
      all_of(set_cols),
      ~ as.logical(as.numeric(.x))
    )
  )

# ============================================================
# FILTRO OPCIONAL DE GERAÇÃO
# ============================================================

if ("Generation" %in% names(upset_input_xlsx) && !identical(generation_to_use, "all")) {
  upset_input_xlsx <- upset_input_xlsx %>%
    filter(Generation == generation_to_use)
}

if (nrow(upset_input_xlsx) == 0) {
  stop("A base ficou vazia após o filtro aplicado.")
}

# ============================================================
# CONVERTENDO PARA FORMATO LONGO
# ============================================================

upset_long <- upset_input_xlsx %>%
  pivot_longer(
    cols = all_of(set_cols),
    names_to = "Model",
    values_to = "Selected"
  ) %>%
  mutate(
    Model = factor(Model, levels = model_order),
    ScenarioLabel = paste0("F = ", format(TargetInbreeding, nsmall = 2))
  )

# ============================================================
# CONTANDO NÚMERO DE INDIVÍDUOS POR MODELO EM CADA UNIDADE
# ============================================================
# Se generation_to_use = "all":
#   conta em cada Replication x Generation x TargetInbreeding
# e depois tira a média desses valores.
#
# Se generation_to_use = 5, por exemplo:
#   conta em cada Replication x TargetInbreeding
# dentro da geração 5 e depois tira a média das repetições.

group_count_cols <- c("TargetInbreeding", "ScenarioLabel", "Replication", "Model")

if ("Generation" %in% names(upset_long) && identical(generation_to_use, "all")) {
  group_count_cols <- c(group_count_cols, "Generation")
}

model_size_by_unit <- upset_long %>%
  group_by(across(all_of(group_count_cols))) %>%
  summarise(
    NumberIndividuals = sum(Selected, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# RESUMO FINAL PARA O GRÁFICO
# ============================================================

model_size_summary <- model_size_by_unit %>%
  group_by(TargetInbreeding, ScenarioLabel, Model) %>%
  summarise(
    MeanNumberIndividuals = mean(NumberIndividuals, na.rm = TRUE),
    SDNumberIndividuals = sd(NumberIndividuals, na.rm = TRUE),
    N = sum(!is.na(NumberIndividuals)),
    SENumberIndividuals = SDNumberIndividuals / sqrt(N),
    ymin = MeanNumberIndividuals - SDNumberIndividuals,
    ymax = MeanNumberIndividuals + SDNumberIndividuals,
    .groups = "drop"
  ) %>%
  mutate(
    Model = factor(Model, levels = model_order)
  )

print(model_size_summary)

# ============================================================
# TEXTO DO SUBTÍTULO
# ============================================================

subtitle_text <- if (identical(generation_to_use, "all")) {
  "Bars represent means across recurrent selection cycles and replications; error bars represent ±1 SD"
} else {
  paste0(
    "Bars represent means across replications for Generation = ",
    generation_to_use,
    "; error bars represent ±1 SD"
  )
}

# ============================================================
# GRÁFICO
# ============================================================

p_model_size_all_f <- ggplot(
  model_size_summary,
  aes(
    x = Model,
    y = MeanNumberIndividuals,
    fill = Model
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
    width = 0.25,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ ScenarioLabel,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_x_discrete(
    labels = pretty_model_name,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean number of selected individuals",
    fill = "Model",
    title = "Mean number of selected individuals by model",
    subtitle = subtitle_text
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
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

print(p_model_size_all_f)

# ============================================================
# SALVANDO O GRÁFICO
# ============================================================

pasta_saida_grafico <- dirname(arquivo_excel)

nome_grafico <- paste0(
  "mean_number_selected_individuals_by_model",
  if (identical(generation_to_use, "all")) {
    "_all_generations"
  } else {
    paste0("_generation_", generation_to_use)
  },
  ".png"
)

arquivo_grafico <- file.path(
  pasta_saida_grafico,
  nome_grafico
)

ggsave(
  filename = arquivo_grafico,
  plot = p_model_size_all_f,
  width = 12,
  height = 5,
  dpi = 500,
  bg = "white"
)

cat("\nGráfico salvo em:\n")
cat(normalizePath(arquivo_grafico), "\n")
