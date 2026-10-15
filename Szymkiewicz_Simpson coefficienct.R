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

f_value <- 0.75

# Use "all" para calcular a média considerando todos os ciclos/gerações.
# Use um número, por exemplo 5, para usar apenas uma geração específica.
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
# LENDO OS DADOS DO XLSX
# ============================================================

containment_base <- openxlsx::read.xlsx(
  xlsxFile = arquivo_excel,
  sheet = "UpSet_Input"
)

# ============================================================
# FUNÇÃO PARA AJUSTAR NOMES NO GRÁFICO
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
# CHECAGENS INICIAIS
# ============================================================

if (!("Replication" %in% names(containment_base))) {
  containment_base <- containment_base %>%
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

missing_required <- setdiff(
  required_cols,
  names(containment_base)
)

if (length(missing_required) > 0) {
  stop(
    paste0(
      "As seguintes colunas obrigatórias não foram encontradas na aba 'UpSet_Input': ",
      paste(missing_required, collapse = ", ")
    )
  )
}

# ============================================================
# FILTRANDO F E, SE QUISER, A GERAÇÃO
# ============================================================

containment_base <- containment_base %>%
  mutate(
    TargetInbreeding = round(as.numeric(TargetInbreeding), 2),
    Replication = as.integer(Replication)
  ) %>%
  filter(
    abs(TargetInbreeding - f_value) < 1e-8
  )

if ("Generation" %in% names(containment_base)) {
  containment_base <- containment_base %>%
    mutate(Generation = as.integer(Generation))
  
  if (!identical(generation_to_use, "all")) {
    containment_base <- containment_base %>%
      filter(Generation == generation_to_use)
  }
}

if (nrow(containment_base) == 0) {
  stop(
    paste0(
      "Não existem linhas na aba 'UpSet_Input' para TargetInbreeding = ",
      format(f_value, nsmall = 2),
      if (!identical(generation_to_use, "all")) {
        paste0(" e Generation = ", generation_to_use)
      } else {
        ""
      },
      "."
    )
  )
}

# ============================================================
# IDENTIFICANDO COLUNAS DOS MODELOS
# ============================================================

set_cols <- setdiff(
  names(containment_base),
  metadata_cols
)

if (length(set_cols) < 2) {
  stop("Menos de dois modelos foram encontrados na aba 'UpSet_Input'.")
}

containment_base <- containment_base %>%
  mutate(
    across(
      all_of(set_cols),
      ~ as.logical(as.numeric(.x))
    )
  )

model_order <- set_cols

# ============================================================
# CRIANDO PARES DE COMPARAÇÃO
# ============================================================

model_pairs <- expand.grid(
  SourceModel = model_order,
  ReferenceModel = model_order,
  stringsAsFactors = FALSE
) %>%
  mutate(
    idxSource = match(SourceModel, model_order),
    idxReference = match(ReferenceModel, model_order)
  )

calc_directional_containment_one_group <- function(data_group, source_model, reference_model) {
  ids_source <- data_group %>%
    filter(.data[[source_model]]) %>%
    pull(IndividualID) %>%
    unique()
  
  ids_reference <- data_group %>%
    filter(.data[[reference_model]]) %>%
    pull(IndividualID) %>%
    unique()
  
  inter_n <- length(intersect(ids_source, ids_reference))
  source_size <- length(ids_source)
  reference_size <- length(ids_reference)
  
  if (source_size == 0) {
    return(
      data.frame(
        DirectionalContainment = NA_real_,
        IntersectionSize = inter_n,
        SourceSetSize = source_size,
        ReferenceSetSize = reference_size
      )
    )
  }
  
  data.frame(
    DirectionalContainment = inter_n / source_size,
    IntersectionSize = inter_n,
    SourceSetSize = source_size,
    ReferenceSetSize = reference_size
  )
}

# ============================================================
# CALCULANDO CONTENÇÃO DIRECIONAL
# ============================================================

group_cols <- c("Replication")

if ("Generation" %in% names(containment_base)) {
  group_cols <- c(group_cols, "Generation")
}

directional_containment <- containment_base %>%
  group_split(across(all_of(group_cols))) %>%
  map_dfr(function(data_group) {
    
    rep_id <- unique(data_group$Replication)[1]
    
    gen_id <- if ("Generation" %in% names(data_group)) {
      unique(data_group$Generation)[1]
    } else {
      NA_integer_
    }
    
    purrr::map_dfr(seq_len(nrow(model_pairs)), function(i) {
      
      source_model <- model_pairs$SourceModel[i]
      reference_model <- model_pairs$ReferenceModel[i]
      
      containment_result <- calc_directional_containment_one_group(
        data_group = data_group,
        source_model = source_model,
        reference_model = reference_model
      )
      
      data.frame(
        SourceModel = source_model,
        ReferenceModel = reference_model,
        idxSource = model_pairs$idxSource[i],
        idxReference = model_pairs$idxReference[i],
        Replication = rep_id,
        Generation = gen_id,
        DirectionalContainment = containment_result$DirectionalContainment,
        IntersectionSize = containment_result$IntersectionSize,
        SourceSetSize = containment_result$SourceSetSize,
        ReferenceSetSize = containment_result$ReferenceSetSize
      )
    })
  })

directional_containment_summary <- directional_containment %>%
  group_by(
    SourceModel,
    ReferenceModel,
    idxSource,
    idxReference
  ) %>%
  summarise(
    MeanContainment = mean(DirectionalContainment, na.rm = TRUE),
    SDContainment = sd(DirectionalContainment, na.rm = TRUE),
    NContainment = sum(!is.na(DirectionalContainment)),
    SEContainment = SDContainment / sqrt(NContainment),
    MeanIntersectionSize = mean(IntersectionSize, na.rm = TRUE),
    MeanSourceSetSize = mean(SourceSetSize, na.rm = TRUE),
    MeanReferenceSetSize = mean(ReferenceSetSize, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    RowModel = factor(SourceModel, levels = rev(model_order)),
    ColModel = factor(ReferenceModel, levels = model_order),
    LabelContainment = ifelse(
      is.na(MeanContainment),
      "",
      sprintf("%.2f", MeanContainment)
    )
  )

print(directional_containment_summary)

# ============================================================
# GRÁFICO
# ============================================================

titulo_geracao <- if (identical(generation_to_use, "all")) {
  "All generations"
} else {
  paste0("Generation = ", generation_to_use)
}

p_directional_containment <- ggplot(
  directional_containment_summary,
  aes(
    x = ColModel,
    y = RowModel,
    fill = MeanContainment
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = LabelContainment),
    size = 3
  ) +
  scale_x_discrete(
    labels = pretty_model_name,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = pretty_model_name,
    drop = FALSE
  ) +
  scale_fill_gradient(
    low = "#F8FAFC",
    high = "#1E293B",
    limits = c(0, 1),
    name = "Mean containment"
  ) +
  coord_equal() +
  labs(
    x = "Reference model",
    y = "Source model",
    title = paste0(
      "Directional containment of selected individuals - F = ",
      format(f_value, nsmall = 2)
    ),
    subtitle = paste0(
      "Data source: Excel sheet 'UpSet_Input' | ",
      titulo_geracao,
      " | Cell value = proportion of individuals from the row model also present in the column model"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

print(p_directional_containment)

# ============================================================
# SALVANDO O GRÁFICO
# ============================================================

pasta_saida_grafico <- dirname(arquivo_excel)

nome_grafico <- paste0(
  "directional_containment_F_",
  gsub("\\.", "_", format(f_value, nsmall = 2)),
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
  plot = p_directional_containment,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)

cat("\nGráfico salvo em:\n")
cat(normalizePath(arquivo_grafico), "\n")
