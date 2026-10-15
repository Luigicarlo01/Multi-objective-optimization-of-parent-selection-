packages <- c(
  "openxlsx",
  "dplyr",
  "ggplot2"
)

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(char = packages)

# ============================================================
# CONFIGURAÇÕES
# ============================================================

pasta_saida_base <- "Resultados_cenarios_endogamia_um_ciclo_sem_repeticoes"

arquivo_excel <- file.path(
  "Resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes",
  "Todos_os_modelos",
  "resultados_cenarios_endogamia_5_ciclos_candidatos300_pairwise_reduzido5x_sem_repeticoes.xlsx"
)

aba_dados <- "Combined_Results"

if (!file.exists(arquivo_excel)) {
  stop(
    paste0(
      "O arquivo Excel não foi encontrado em:\n",
      arquivo_excel,
      "\n\nConfira se o caminho em 'arquivo_excel' está correto."
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
# LENDO OS DADOS DIRETAMENTE DO XLSX
# ============================================================

dados <- openxlsx::read.xlsx(
  xlsxFile = arquivo_excel,
  sheet = aba_dados
)

# ============================================================
# FUNÇÃO PARA PADRONIZAR NOMES DOS ALGORITMOS
# ============================================================

clean_algorithm_name <- function(x) {
  x <- as.character(x)
  
  dplyr::case_when(
    grepl("^Strategy[_ ]1", x) ~ "Strategy 1",
    grepl("^Strategy[_ ]2", x) ~ "Strategy 2",
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

algorithm_levels <- c(
  "Strategy 1",
  "Strategy 2",
  "EA",
  "Local",
  "SA",
  "GA",
  "GRASP",
  "TABU",
  "EDA",
  "PSO"
)

# ============================================================
# CHECANDO COLUNAS NECESSÁRIAS
# ============================================================

required_cols <- c(
  "TargetInbreeding",
  "Strategy",
  "GenomicEndogamyF",
  "MeanPhenotypicValue"
)

missing_cols <- setdiff(required_cols, names(dados))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "As seguintes colunas obrigatórias não foram encontradas na aba '",
      aba_dados,
      "': ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# PREPARANDO OS DADOS
# ============================================================

dados_plot <- dados %>%
  mutate(
    TargetInbreeding = round(as.numeric(TargetInbreeding), 2),
    GenomicEndogamyF = as.numeric(GenomicEndogamyF),
    MeanPhenotypicValue = as.numeric(MeanPhenotypicValue),
    Algorithm = clean_algorithm_name(Strategy),
    Algorithm = factor(Algorithm, levels = algorithm_levels),
    ScenarioLabel = paste0(
      "Initial F = ",
      sprintf("%.2f", TargetInbreeding)
    )
  ) %>%
  filter(
    TargetInbreeding != 0
  )

if (nrow(dados_plot) == 0) {
  stop("A base ficou vazia após remover TargetInbreeding = 0.")
}

scenario_levels <- dados_plot %>%
  distinct(TargetInbreeding, ScenarioLabel) %>%
  arrange(TargetInbreeding) %>%
  pull(ScenarioLabel)

dados_plot <- dados_plot %>%
  mutate(
    ScenarioLabel = factor(
      ScenarioLabel,
      levels = scenario_levels
    )
  )

# ============================================================
# RESUMO PARA O GRÁFICO
# ============================================================

resumo_xy <- dados_plot %>%
  group_by(
    ScenarioLabel,
    TargetInbreeding,
    Algorithm
  ) %>%
  summarise(
    MeanEndogamy = mean(GenomicEndogamyF, na.rm = TRUE),
    SDEndogamy = sd(GenomicEndogamyF, na.rm = TRUE),
    MeanPheno = mean(MeanPhenotypicValue, na.rm = TRUE),
    SDPheno = sd(MeanPhenotypicValue, na.rm = TRUE),
    NEndogamy = sum(!is.na(GenomicEndogamyF)),
    NPheno = sum(!is.na(MeanPhenotypicValue)),
    .groups = "drop"
  ) %>%
  mutate(
    SDEndogamy = ifelse(is.na(SDEndogamy), 0, SDEndogamy),
    SDPheno = ifelse(is.na(SDPheno), 0, SDPheno),
    xmin = MeanEndogamy - SDEndogamy,
    xmax = MeanEndogamy + SDEndogamy,
    ymin = MeanPheno - SDPheno,
    ymax = MeanPheno + SDPheno
  )

print(resumo_xy)

# ============================================================
# AJUSTE DOS TAMANHOS DAS BARRAS DE ERRO
# ============================================================

x_range <- diff(range(resumo_xy$MeanEndogamy, na.rm = TRUE))
y_range <- diff(range(resumo_xy$MeanPheno, na.rm = TRUE))

if (!is.finite(x_range) || x_range == 0) x_range <- 1
if (!is.finite(y_range) || y_range == 0) y_range <- 1

cap_x <- 0.015 * x_range
cap_y <- 0.015 * y_range

# ============================================================
# GRÁFICO
# ============================================================

p_cross <- ggplot(
  resumo_xy,
  aes(
    x = MeanEndogamy,
    y = MeanPheno,
    color = Algorithm
  )
) +
  geom_segment(
    aes(
      x = xmin,
      xend = xmax,
      y = MeanPheno,
      yend = MeanPheno
    ),
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_segment(
    aes(
      x = MeanEndogamy,
      xend = MeanEndogamy,
      y = ymin,
      yend = ymax
    ),
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_segment(
    aes(
      x = xmin,
      xend = xmin,
      y = MeanPheno - cap_y,
      yend = MeanPheno + cap_y
    ),
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_segment(
    aes(
      x = xmax,
      xend = xmax,
      y = MeanPheno - cap_y,
      yend = MeanPheno + cap_y
    ),
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_segment(
    aes(
      x = MeanEndogamy - cap_x,
      xend = MeanEndogamy + cap_x,
      y = ymin,
      yend = ymin
    ),
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_segment(
    aes(
      x = MeanEndogamy - cap_x,
      xend = MeanEndogamy + cap_x,
      y = ymax,
      yend = ymax
    ),
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ ScenarioLabel,
    scales = "free"
  ) +
  labs(
    x = "Genomic endogamy",
    y = "Mean phenotypic value",
    color = "Algorithm",
    title = "Phenotypic value versus genomic endogamy by algorithm",
    subtitle = "Each point is the mean of the algorithm; horizontal and vertical arms represent ±1 SD"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_cross)

# ============================================================
# SALVANDO O GRÁFICO
# ============================================================

pasta_saida_grafico <- dirname(arquivo_excel)

arquivo_grafico <- file.path(
  pasta_saida_grafico,
  "phenotypic_value_vs_genomic_endogamy_by_algorithm.png"
)

ggsave(
  filename = arquivo_grafico,
  plot = p_cross,
  width = 10,
  height = 7,
  dpi = 500,
  bg = "white"
)

cat("\nGráfico salvo em:\n")
cat(normalizePath(arquivo_grafico), "\n")
