library("tidyverse")
library("edgeR")
library("readxl")
library("modelr")
library("wesanderson")
library("patchwork")
library("data.table")
library("eulerr")
library("gtools")
library("rcompanion")

##########################
## BRAF Mutation Analysis
##########################

## Variables
## ----------

comparisons <- tibble(
  gene_1_name = c("NEUROG3", "CHGA", "MUC2", "LGR5"),
  gene_2_name = c("INSM1", "SST", "FCGBP", "AXIN2"),
  gene_1_id = c(
    "ENSG00000122859", "ENSG00000100604",
    "ENSG00000198788", "ENSG00000139292"
  ),
  gene_2 = c(
    "ENSG00000173404", "ENSG00000157005",
    "ENSG00000275395", "ENSG00000168646"
  )
)

## Prepare Mutation Data
## ----------

## Get file and corresponding gene name.

mut_files <- list.files("references", pattern = ".*\\.tsv$", full.names = TRUE) %>%
  keep(~ !str_detect(., "Mution"))

sample_names <- mut_files %>%
  str_extract("(?<=/)[[:alnum:]]+(?=\\.|-)") %>%
  str_to_upper

names(mut_files) <- sample_names

## Clean and combine data.

mut_data <- map(mut_files, function(dataset) {
  muts <- dataset %>%
    read.delim(header = FALSE, skip = 1, stringsAsFactors = FALSE) %>%
    as_tibble

  if (ncol(muts) == 24) {
    muts <- bind_rows(select(muts, 1:12), select(muts, 13:24))
  } else if (ncol(muts) == 36) {
    muts <- bind_rows(select(muts, 1:12), select(muts, 13:24), select(muts ,25:36))
  }
  muts <- select(muts, seq_len(12))

  colnames(muts) <- c(
    "sample", "chr", "start", "end", "gene", "reference",
    "alt", "altGene", "effect", "aminoAcid", "rnaVaf", "dnaVaf"
  )
  muts <- filter_all(muts, all_vars(!is.na(.)))

  return(muts)
}) %>%
  bind_rows(.id = "dataset") %>%
  mutate(sample = str_replace_all(sample, "-", "."))

## Keep only useful mutations.

mut_data <- filter(
  mut_data,
  !(effect %in% c("no variant", "synonymous_variant")),
  gene != "" & !is.na(gene)
)

fwrite(
  mut_data, file.path("audit", "mut_data_raw.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

mut_data %>%
  count(dataset, gene, name="sample_count") %>%
  fwrite(
    file.path("audit", "mut_data_counts.tsv"), sep="\t",
    col.names=TRUE, row.names=FALSE, quote=FALSE
  )

## Prepare Expression Data
## -----------

exp_file <- file.path("references", "expression", "TCGA-COAD.htseq_counts.tsv")

## Load counts data.

counts <- exp_file %>%
  read.delim(header = TRUE, sep = "\t", stringsAsFactor = FALSE) %>%
  as_tibble(.name_repair = "unique") %>%
  filter(!grepl(Ensembl_ID, pattern = "^__")) %>%
  mutate(Ensembl_ID = str_replace(Ensembl_ID, "\\.[0-9]+$", "")) %>%
  mutate_if(is.numeric, ~ {(2^.) + 1})

## Counts were log2(x+1) transformed, so revert them back and TMM normalize the counts.

norm_counts <- counts %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix %>%
  DGEList %>%
  calcNormFactors %>%
  cpm %>%
  as.data.frame %>%
  as_tibble(.name_repair = "unique", rownames = "Ensembl_ID")

quants <- norm_counts %>%
  pivot_longer(starts_with("TCGA"), names_to="sample", values_to="exp") %>%
  pivot_wider(names_from="Ensembl_ID", values_from="exp") %>%
  transmute(
    sample = sample,
    NEUROG3_quants = factor(ntile(ENSG00000122859, 5)),
    INSM1_quants = factor(ntile(ENSG00000173404, 5)),
    SST_quants = factor(ntile(ENSG00000157005, 5)),
    CHGA_quants = factor(ntile(ENSG00000100604, 5)),
    FCGBP_quants = factor(ntile(ENSG00000275395, 5)),
    MUC2_quants = factor(ntile(ENSG00000198788, 5)),
    AXIN2_quants = factor(ntile(ENSG00000168646, 5)),
    LGR5_quants = factor(ntile(ENSG00000139292, 5))
  )

## Merge the mutation data back in.

muts_to_merge <- mut_data %>%
  select(sample, gene) %>%
  distinct

merged <- quants %>%
  left_join(muts_to_merge, by = "sample") %>%
  filter(gene %in% c("APC", "BRAF", "KRAS", "PIK3CA", "SMAD4", "TP53")) %>%
  split(., .$gene)

## Permutation Test
## ----------

ncounts <- as.data.table(norm_counts)
ncounts <- melt(ncounts, id.vars="Ensembl_ID", variable.name="sample", value.name="exp")
ncounts <- dcast(ncounts, sample ~ Ensembl_ID, value.var="exp")
ncounts <- ncounts[, .(sample,
  NEUROG3_exp=ENSG00000122859, INSM1_exp=ENSG00000173404,
  SST_exp=ENSG00000157005, CHGA_exp=ENSG00000100604,
  FCGBP_exp=ENSG00000275395, MUC2_exp=ENSG00000198788,
  AXIN2_exp=ENSG00000168646, LGR5_exp=ENSG00000139292
)]

nmuts <- as.data.table(mut_data)[, .(sample, gene)]
nmuts <- unique(nmuts)
nmuts <- nmuts[gene %in% c("APC", "BRAF", "KRAS", "PIK3CA", "SMAD4", "TP53")]

nmerged <- merge(ncounts, nmuts, by="sample")

fwrite(
  nmerged, file.path("audit", "mut_data_expression.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

cols <- str_subset(colnames(nmerged), "_exp")
nmerged[, (cols) := lapply(.SD, ntile, 5), .SDcols=cols]

fwrite(
  nmerged, file.path("audit", "mut_data_quintiles.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Observed

obs_results <- pmap(comparisons, function(...) {

  args <- list(...)
  gene_1 <- str_c(args$gene_1_name, "_exp")
  gene_2 <- str_c(args$gene_2_name, "_exp")
  genes <- c(gene_1, gene_2)

  cols <- c("gene", genes, "sample")
  observed <- unique(nmerged[, ..cols])
  observed[, sample := NULL]
  setnames(observed, old=genes, new=c("gene_1", "gene_2"))
  n_observed <- observed[, .(n_obs = .N), by = .(gene_1, gene_2, gene)]
  n_observed <- n_observed[
    CJ(gene_1 = seq_len(5), gene_2 = seq_len(5), gene = unique(n_observed$gene)),
    on = .(gene_1, gene_2, gene)
  ]
  setnafill(n_observed, cols="n_obs", fill=0)
  setnames(n_observed, old=c("gene_1", "gene_2"), new=genes)
  return(n_observed)
})

names(obs_results) <- comparisons %>%
  transmute(comparison=str_c(gene_1_name, "vs", gene_2_name, sep="_")) %>%
  pull(comparison)

iwalk(
  obs_results,
  ~fwrite(
    .x, file.path("audit", str_c(.y, ".tsv")), sep="\t",
    col.names=TRUE, row.names=FALSE, quote=FALSE
  )
)

obs_results %>%
  map(~group_by(., gene) %>% summarize(count=sum(n_obs))) %>%
  bind_rows(.id="gene_pair") %>%
  fwrite(
    file.path("audit", "gene_pair_totals.tsv"), sep="\t",
    col.names=TRUE, row.names=FALSE, quote=FALSE
  )

wes_colors <- as.character(wes_palette("Zissou1", 100, "continuous"))

plots <- map(obs_results, function(x) {
  comparisons <- str_subset(colnames(x), "_exp$")
  if (any(str_detect(comparisons, "NEUROG3"))) {
    comparisons <- rev(comparisons)
  }

  p <- ggplot(x, aes_string(x=comparisons[1], y=comparisons[2])) +
    geom_tile(aes(fill=n_obs, color=n_obs)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text=element_text(size=18)
    ) +
    scale_fill_gradientn(colors = wes_colors) +
    scale_color_gradientn(colors = wes_colors) +
    coord_fixed() +
    facet_wrap(gene ~ ., ncol=3)
  return(p)
})

plots <- wrap_plots(plots)

ggsave(
  file.path("results", "mutation_tables", "mutation_plot_observed.pdf"),
  plot=plots, height=12, width=12, device=cairo_pdf
)

## Simulated

n <- 10000
perm_results <- pmap(comparisons, function(...) {

  args <- list(...)
  gene_1 <- str_c(args$gene_1_name, "_exp")
  gene_2 <- str_c(args$gene_2_name, "_exp")
  genes <- c(gene_1, gene_2)

  cols <- c("gene", genes, "sample")
  observed <- unique(nmerged[, ..cols])
  observed <- nmerged[, ..cols]
  observed[, sample := NULL]
  setnames(observed, old=genes, new=c("gene_1", "gene_2"))
  n_observed <- observed[, .(n_obs = .N), by = .(gene_1, gene_2, gene)]
  n_observed <- n_observed[
    CJ(gene_1 = seq_len(5), gene_2 = seq_len(5), gene = unique(n_observed$gene)),
    on = .(gene_1, gene_2, gene)
  ]
  setnafill(n_observed, cols="n_obs", fill=0)

  permuted <- lapply(seq_len(n), function(x) {
    perm <- observed[, .(gene_1 = sample(gene_1), gene_2 = sample(gene_2), gene)]
    perm <- perm[, .(n_sim = .N), by = .(gene_1, gene_2, gene)]
    return(perm)
  })
  permuted <- rbindlist(permuted)
  permuted <- merge(permuted, n_observed, by=c("gene_1", "gene_2", "gene"), all=TRUE)
  setnafill(permuted, cols=c("n_sim"), fill=0)
  permuted[, c("pval", "diff") := list(
    n_sim >= n_obs,
    log2(n_obs + 1) - log2(n_sim + 1)
  )]
  permuted <- permuted[,
    .(pval = (sum(pval) + 1) / (n + 1), mean_diff = mean(diff)),
    by = .(gene_1, gene_2, gene)
  ]
  permuted[, FDR := p.adjust(pval, "fdr")]
  setnames(permuted, old=c("gene_1", "gene_2"), new=genes)
  return(permuted)
})

names(perm_results) <- comparisons %>%
  transmute(comparison=str_c(gene_1_name, "vs", gene_2_name, sep="_")) %>%
  pull(comparison)

if (!dir.exists(file.path("results", "mutation_tables"))) {
  dir.create(file.path("results", "mutation_tables"), recursive=TRUE)
}

iwalk(perm_results, function(x, y) {
  fwrite(
    x, file.path("results", "mutation_tables", str_c(y, ".tsv")),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )
})

saveRDS(perm_results, file.path("results", "mutation_tables", "perm_results.RDS"))

## Plot the results.

perm_results <- map(perm_results, function(x) {
  x[, sig := ifelse(FDR < 0.05, "FDR < 0.05", "n.s.")]
  x[, sig := factor(sig, levels=c("FDR < 0.05", "n.s."))]
  cols <- str_subset(colnames(x), "_exp$")
  x[, (cols) := lapply(.SD, function(x) factor(x, levels=seq_len(5))), .SDcols=cols]
  x[, mean_diff_sig := ifelse(FDR < 0.05, "*", "")]
  return(x)
})


# Significance plot.
wes_colors <- as.character(wes_palette("Zissou1", 2, "continuous"))

plots <- map(perm_results, function(x) {
  comparisons <- str_subset(colnames(x), "_exp$")
  if (any(str_detect(comparisons, "NEUROG3"))) {
    comparisons <- rev(comparisons)
  }

  p <- ggplot(x, aes_string(x=comparisons[1], y=comparisons[2])) +
    geom_tile(aes(fill=sig), color="white", width=0.95, height=0.95) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text=element_text(size=18)
    ) +
    scale_fill_manual(values = c("FDR < 0.05"=wes_colors[2], "n.s."="grey")) +
    coord_fixed() +
    facet_wrap(gene ~ ., ncol=3)
  return(p)
})

plots <- wrap_plots(plots)

ggsave(
  file.path("results", "mutation_tables", "mutation_sig_plot.pdf"),
  plot=plots, height=12, width=12, device=cairo_pdf
)

# FDR and mean log2 difference plot.
wes_colors <- as.character(wes_palette("Zissou1", 100, "continuous"))

plots <- map(perm_results, function(x) {
  comparisons <- str_subset(colnames(x), "_exp$")
  if (any(str_detect(comparisons, "NEUROG3"))) {
    comparisons <- rev(comparisons)
  }

  p <- ggplot(x, aes_string(x=comparisons[1], y=comparisons[2])) +
    geom_tile(aes(fill=mean_diff, color=mean_diff)) +#, color="white", width=0.95, height=0.95) +
    geom_text(aes(label=mean_diff_sig), nudge_y=-0.2, color="white", size=8) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text=element_text(size=18)
    ) +
    scale_fill_gradientn(colors = wes_colors) +
    scale_color_gradientn(colors = wes_colors) +
    coord_fixed() +
    facet_wrap(gene ~ ., ncol=3)
  return(p)
})

plots <- wrap_plots(plots)

ggsave(
  file.path("results", "mutation_tables", "mutation_FDR_plot.pdf"),
  plot=plots, height=12, width=12, device=cairo_pdf
)


################
## BRAF V600E ##
################

braf_muts <- mut_data %>%
  select(sample, gene, amino_acid=aminoAcid) %>%
  filter(gene %in% c("BRAF", "KRAS", "PIk3CA")) %>%
  left_join(quants, by="sample") %>%
  mutate(V600E=case_when(
    gene != "BRAF" ~ "not_BRAF",
    gene == "BRAF" & str_detect(amino_acid, "V600E") ~ "V600E",
    gene == "BRAF" & !str_detect(amino_acid, "V600E") ~ "other"
  ))

###########################
## BRAF Mutation Domains ##
###########################

## Load and prepare mutation data.

braf_domain_regions <- read_xlsx(file.path("references", "BRAF_domains_Dankner_Oncogene_2018.xlsx"))
conserved_regions <- filter(braf_domain_regions, Domain %in% c("CR1", "CR2", "CR3_Kinase"))
defined_domains <- filter(braf_domain_regions, !Domain %in% c("CR1", "CR2", "CR3_Kinase"))

## Add mutation data back to braf_muts.

# Get the mutated amino-acid.
braf_domains <- braf_muts %>%
  pull(amino_acid) %>%
  str_extract_all("\\d+", simplify = TRUE) %>%
  as_tibble %>%
  dplyr::rename("start" = 1, "end" = 2) %>%
  bind_cols(braf_muts, .) %>%
  mutate(end = ifelse(end == "", start, end))

# Define the domain the mutation is in.
braf_domains <- braf_domains %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  mutate(conserved_regions = case_when(
    gene != "BRAF" ~ "not_BRAF",
    start >= 150 & start <= 290 ~ "CR1",
    start >= 360 & start <= 375 ~ "CR2",
    start >= 457 & start <= 717 ~ "CR3_Kinase",
    TRUE ~ "other"
  )) %>%
  mutate(defined_domains = case_when(
    gene != "BRAF" ~ "not_BRAF",
    start >= 155 & start <= 227 ~ "RBD",
    start >= 464 & start <= 471 ~ "P_loop",
    start >= 492 & start < 504 ~ "aC_helix",
    start >= 504 & start <= 511 ~ "DIF",
    start >= 574 & start <= 581 ~ "CL",
    start >= 594 & start <= 596 ~ "DFG",
    start > 596 & start <= 623 ~ "AS",
    TRUE ~ "other"
  ))

##############
## Subtypes ##
##############

## Normal vs. Cancer.

norm_split <- read.delim(
    file.path("references", "Mution_normal_annotation.tsv"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  ) %>%
  dplyr::rename("sample_type" = 2) %>%
  mutate(
    sample = str_replace_all(sample, "-", "."),
    sample_type = str_replace_all(sample_type, " ", "_")
  )

subtype_merged <- left_join(braf_domains, norm_split, by="sample")

## Filter out data not being used for this analysis.

subtype <- subtype_merged %>%
  filter(conserved_regions != "other") %>%
  filter(across(ends_with("quants"), ~!is.na(.x)))

## Define Subtype.

subtype <- subtype %>%
  group_by(sample) %>%
  mutate(
    KRAS_EP = ifelse(
      NEUROG3_quants %in% c(4,5) & INSM1_quants %in% c(4,5) & any(gene == "KRAS"),
      "*", ""
    ),
    KRAS_EC = ifelse(
      SST_quants %in% c(4,5) & CHGA_quants %in% c(4,5) & any(gene == "KRAS"),
      "*", ""
    ),
    BRAF_EP = ifelse(
      NEUROG3_quants %in% c(4,5) & INSM1_quants %in% c(4,5) & any(gene == "BRAF"),
      "*", ""
    ),
    BRAF_EC = ifelse(
      SST_quants %in% c(4,5) & CHGA_quants %in% c(4,5) & any(gene == "BRAF"),
      "*", ""
    )
  ) %>%
  ungroup

dir.create(file.path("audit", "subtype"))

fwrite(
  subtype, file.path("audit", "subtype", "subtypes.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Euler plot.

euler_plot <- subtype %>%
  select(matches("(EP|EC)$")) %>%
  mutate(across(everything(), ~.x == "*")) %>%
  euler %>%
  plot(quantities=TRUE, fill=as.character(wes_palette("Zissou1", 4, "continuous")))

cairo_pdf(file.path("audit", "subtype", "subtype_euler.pdf"), height=4, width=4)
euler_plot; dev.off()

###################
## Clinical Data ##
###################

## Prepare sample type.

clinical_type <- norm_split

## Prepare subtype data.

clinical_subtypes <- subtype %>%
  rename(mut_gene=gene) %>%
  select(sample, matches("(EP|EC)$")) %>%
  distinct

## Preparing count data.

clinical_counts <- norm_counts %>%
  filter(Ensembl_ID %in% c("ENSG00000172238", "ENSG00000198788")) %>%
  pivot_longer(!Ensembl_ID, names_to="sample", values_to="exp") %>%
  rename(gene=Ensembl_ID) %>%
  mutate(gene=fct_recode(gene, ATOH1="ENSG00000172238", MUC2="ENSG00000198788"))

## Merging the clinical type, subtype, and counts data.

clinical <- reduce(list(clinical_type, clinical_subtypes, clinical_counts), full_join, by="sample") %>%
  as_tibble %>%
  filter(across(c(gene, exp), ~!is.na(.x))) %>%
  mutate(sample_simple=str_replace(sample, "\\.[[:alnum:]]{3}$", "")) %>%
  group_by(normal_tissue=str_detect(sample_type, "Normal")) %>%
  nest
  
## Importing clinical data.

clinical_data <- file.path("references", "clinical", "clinical.tsv") %>%
  read_tsv %>%
  select(primary_diagnosis, sample=submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  distinct

## Merge clinical data into subtype data.

clinical <- clinical %>%
  rowwise %>%
  mutate(data=list(
    if (normal_tissue) {
      data 
    } else {
      left_join(data, clinical_data, by=c("sample_simple"="sample"))
    }
  )) %>%
  unnest(cols=data) %>%
  ungroup %>%
  distinct %>%
  select(!normal_tissue) %>%
  group_by(sample_simple) %>%
  filter(n() <= 4) %>%
  ungroup

## Define subtype.

clinical_sub <- clinical %>%
  rowwise %>%
  mutate(subtype=case_when(
    str_detect(sample_type, "Normal") ~ "normal_matched",
    any(is.na(c_across(matches("(EP|EC)$")))) ~ "other",
    all(c_across(matches("(EP|EC)$")) == "") ~ "other",
    sum(c_across(matches("(EP|EC)$")) == "*") > 1 ~ "mixed_subtype",
    KRAS_EP == "*" ~ "KRAS_EP",
    KRAS_EC == "*" ~ "KRAS_EC",
    BRAF_EP == "*" ~ "BRAF_EP",
    BRAF_EC == "*" ~ "BRAF_EC"
  )) %>%
  ungroup %>%
  mutate(
    exp=log2(exp + 1),
    subtype=factor(subtype, levels=c(
      "normal_matched", "BRAF_EP", "BRAF_EC", "KRAS_EP",
      "KRAS_EC", "mixed_subtype", "other"
    )),
    primary_diagnosis=ifelse(subtype == "normal_matched", "Normal", primary_diagnosis)
  ) %>%
  filter(!is.na(primary_diagnosis))

fwrite(
  clinical_sub, file.path("audit", "subtype", "subtype_clinical.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Plot.

p <- ggplot(clinical_sub, aes(x=subtype, y=exp)) +
  geom_jitter(aes(color=primary_diagnosis), width=0.2, size=0.5) +
  geom_boxplot(fill=NA, width=0.5, outlier.shape=NA) +
  facet_wrap(gene~., scales="free_y", ncol=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(
  file.path("audit", "subtype", "subtype_clinical_plot.pdf"),
  plot=p, device=cairo_pdf, height=6, width=8
)

###########################
## Subtype Pairwise Test ##
###########################

## Pairwise all samples.

pvals_inter <- clinical_sub %>%
  split(., .$gene) %>%
  map(function(x) {
    x %>%
      {pairwise.wilcox.test(.$exp, .$subtype, p.adjust="BH")} %>%
      pluck("p.value") %>%
      as_tibble(rownames="set_1") %>%
      pivot_longer(!set_1, names_to="set_2", values_to="padj") %>%
      filter(!is.na(padj))
  }) %>%
  bind_rows(.id="gene")

fwrite(
  pvals_inter, file.path("audit", "subtype", "subtype_pairwise_all.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Pairwise conditions between samples.

pvals_intra <- clinical_sub %>%
  group_by(gene, subtype) %>%
  nest %>%
  filter(subtype != "normal_matched") %>%
  rowwise %>%
  summarize(pval=list(data %>%
    {pairwise.wilcox.test(.$exp, .$primary_diagnosis, p.adjust="BH")} %>%
    pluck("p.value") %>%
    as_tibble(rownames="set_1") %>%
    pivot_longer(!set_1, names_to="set_2", values_to="padj") %>%
    filter(!is.na(padj))
  )) %>%
  unnest(cols=pval)

fwrite(
  pvals_intra, file.path("audit", "subtype", "subtype_pairwise_within.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## Number in each category.

count_table <- clinical_sub %>%
  select(sample, primary_diagnosis, subtype) %>%
  distinct %>%
  group_by(subtype) %>%
  add_count(name="n_samples_per_subtype") %>%
  group_by(subtype, primary_diagnosis) %>%
  add_count(name="n_diagnosis_per_subtype") %>%
  group_by(primary_diagnosis) %>%
  add_count(name="n_diagnosis_total") %>%
  ungroup %>%
  select(!sample) %>%
  distinct %>%
  complete(subtype, primary_diagnosis, fill=list(n_diagnosis_per_subtype=0)) %>%
  group_by(subtype) %>%
  mutate(n_samples_per_subtype=max(n_samples_per_subtype, na.rm=TRUE)) %>%
  group_by(primary_diagnosis) %>%
  mutate(n_diagnosis_total=max(n_diagnosis_total, na.rm=TRUE)) %>%
  ungroup

count_table %>%
  select(
    subtype, n_samples_per_subtype,
    primary_diagnosis, n_diagnosis_per_subtype,
    n_diagnosis_total
  ) %>%
  arrange(subtype, primary_diagnosis) %>%
  fwrite(
    file.path("audit", "subtype", "diagnosis_count_table.tsv"),
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"
  )

pvals <- count_table %>%
  group_by(primary_diagnosis) %>%
  nest %>%
  rowwise %>%
  summarize(pvals=list(
    data %>%
      mutate(not_diagnosis = n_samples_per_subtype - n_diagnosis_per_subtype) %>%
      select(subtype, n_diagnosis_per_subtype, not_diagnosis) %>%
      column_to_rownames("subtype") %>%
      as.matrix %>%
      pairwiseNominalIndependence(fisher=TRUE, gtest=FALSE, chisq=FALSE, digits=6)
  )) %>%
  unnest(cols=pvals) %>%
  separate(Comparison, into=c("subtype_1", "subtype_2"), sep=" : ")

fwrite(
  pvals, file.path("audit", "subtype", "diagnosis_count_pvals.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

#########################
## EEC Marker Analysis ##
#########################

## Prepare markers to analyze.

genes <- c(
        "ENSG00000173404" = "INSM1",
        "ENSG00000122859" = "NEUROG3",
        "ENSG00000162992" = "NEUROD1",
        "ENSG00000004848" = "ARX",
        "ENSG00000162761" = "LMX1A",
        "ENSG00000106331" = "PAX4"
)

eec_marker_file <- read_xlsx(file.path("references", "EEC_cell_makers_ensmbl_IDs.xlsx"), col_names = FALSE)
eec_markers <- pull(eec_marker_file, 1)
names(eec_markers) <- pull(eec_marker_file, 2)
eec_markers <- eec_markers[eec_markers %in% c("CHGA", "CHGB", "CPE", "PYY", "TAC1")]

genes <- c(genes, eec_markers)

## Prepare data for analysis.

eec_marker_exp <- norm_counts %>%
        gather(key = "sample", value = "exp", -Ensembl_ID) %>%
        filter(Ensembl_ID %in% names(genes)) %>%
        rename(gene = Ensembl_ID) %>%
        select(sample, gene, exp) %>%
        mutate(gene = genes[gene])

eec_marker_subtype <- left_join(subtype_clinical, eec_marker_exp, by = "sample") %>%
        filter(subtype_id %in% c("BRAF_EC", "BRAF_EP", "KRAS_EC", "KRAS_EP", "Normal", "Other")) %>%
        mutate(subtype_id = factor(
                subtype_id,
                levels = c("Normal", "BRAF_EP", "BRAF_EC", "KRAS_EP", "KRAS_EC", "Other")
        )) %>%
        select(sample, subtype_id, gene, exp)

## Pairwise wilcoxon tests.

marker_subtype_stats <- eec_marker_subtype %>%
  split(., .$gene) %>%
  map(function(x) {
    res <- pairwise.wilcox.test(x$exp, x$subtype_id, p.adjust="BH")
    res <- res$p.value %>%
      as_tibble(rownames="subtype_1") %>%
      pivot_longer(-subtype_1, names_to="subtype_2", values_to="padj") %>%
      drop_na %>%
      mutate(sig = ifelse(padj < 0.05, "*", ""))
    return(res)
  }) %>%
  bind_rows(.id="gene")

fwrite(
  marker_subtype_stats, file.path("results", "marker_subtype_stats.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
