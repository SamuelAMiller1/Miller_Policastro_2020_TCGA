library("tidyverse")
library("edgeR")
library("readxl")
library("wesanderson")
library("patchwork")
library("data.table")

##########################
## BRAF Mutation Analysis
##########################

## Variables
## ----------

comparisons <- tibble(
  gene_1_name = c("DCLK1", "DEFA5"),
  gene_2_name = c("PTGS1", "DEFA6"),
  gene_1_id = c("ENSG00000133083", "ENSG00000164816"),
  gene_2 = c("ENSG00000095303", "ENSG00000164822")
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

#mut_data %>%
#  count(dataset, gene, name="sample_count") %>%
#  fwrite(
#    file.path("audit", "mut_data_counts.tsv"), sep="\t",
#    col.names=TRUE, row.names=FALSE, quote=FALSE
#  )

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
    DCLK1_quants = factor(ntile(ENSG00000133083, 5)),
    PTGS1_quants = factor(ntile(ENSG00000095303, 5)),
    DEFA5_quants = factor(ntile(ENSG00000164816, 5)),
    DEFA6_quants = factor(ntile(ENSG00000164822, 5))
  )

## Merge the mutation data back in.

muts_to_merge <- mut_data %>%
  select(sample, gene) %>%
  distinct

merged <- quants %>%
  left_join(muts_to_merge, by = "sample") %>%
  filter(gene == "BRAF") %>%
  split(., .$gene)

## Permutation Test
## ----------

ncounts <- as.data.table(norm_counts)
ncounts <- melt(ncounts, id.vars="Ensembl_ID", variable.name="sample", value.name="exp")
ncounts <- dcast(ncounts, sample ~ Ensembl_ID, value.var="exp")
ncounts <- ncounts[, .(sample,
  DCLK1_exp=ENSG00000133083, PTGS1_exp=ENSG00000095303,
  DEFA5_exp=ENSG00000164816, DEFA6_exp=ENSG00000164822
)]

nmuts <- as.data.table(mut_data)[, .(sample, gene)]
nmuts <- unique(nmuts)
nmuts <- nmuts[gene == "BRAF"]

nmerged <- merge(ncounts, nmuts, by="sample")

#fwrite(
#  nmerged, file.path("audit", "mut_data_expression.tsv"),
#  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
#)

cols <- str_subset(colnames(nmerged), "_exp")
nmerged[, (cols) := lapply(.SD, ntile, 5), .SDcols=cols]

#fwrite(
#  nmerged, file.path("audit", "mut_data_quintiles.tsv"),
#  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
#)

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

#iwalk(
#  obs_results,
#  ~fwrite(
#    .x, file.path("audit", str_c(.y, ".tsv")), sep="\t",
#    col.names=TRUE, row.names=FALSE, quote=FALSE
#  )
#)

#obs_results %>%
#  map(~group_by(., gene) %>% summarize(count=sum(n_obs))) %>%
#  bind_rows(.id="gene_pair") %>%
#  fwrite(
#    file.path("audit", "gene_pair_totals.tsv"), sep="\t",
#    col.names=TRUE, row.names=FALSE, quote=FALSE
#  )

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
  file.path("results", "mutation_tables", "mutation_plot_observed_tuft_paneth.pdf"),
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
    x, file.path("results", "mutation_tables", str_c(y, "_tuft_paneth.tsv")),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )
})

saveRDS(perm_results, file.path("results", "mutation_tables", "perm_results_tuft_paneth.RDS"))

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
#wes_colors <- as.character(wes_palette("Zissou1", 2, "continuous"))
#
#plots <- map(perm_results, function(x) {
#  comparisons <- str_subset(colnames(x), "_exp$")
#  if (any(str_detect(comparisons, "NEUROG3"))) {
#    comparisons <- rev(comparisons)
#  }
#
#  p <- ggplot(x, aes_string(x=comparisons[1], y=comparisons[2])) +
#    geom_tile(aes(fill=sig), color="white", width=0.95, height=0.95) +
#    theme_minimal() +
#    theme(
#      panel.grid = element_blank(),
#      text=element_text(size=18)
#    ) +
#    scale_fill_manual(values = c("FDR < 0.05"=wes_colors[2], "n.s."="grey")) +
#    coord_fixed() +
#    facet_wrap(gene ~ ., ncol=3)
#  return(p)
#})
#
#plots <- wrap_plots(plots)
#
#ggsave(
#  file.path("results", "mutation_tables", "mutation_sig_plot.pdf"),
#  plot=plots, height=12, width=12, device=cairo_pdf
#)

# FDR and mean log2 difference plot.
wes_colors <- as.character(wes_palette("Zissou1", 100, "continuous"))

plots <- map(perm_results, function(x) {
  comparisons <- str_subset(colnames(x), "_exp$")

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
  file.path("results", "mutation_tables", "mutation_FDR_plot_tuft_paneth.pdf"),
  plot=plots, height=8, width=8, device=cairo_pdf
)

#####################
## Regression Plot ##
#####################

comparisons <- tibble(
  gene_1_name = c("DCLK1", "DEFA5", "NEUROG3", "CHGA", "MUC2", "LGR5"),
  gene_2_name = c("PTGS1", "DEFA6", "INSM1", "SST", "FCGBP", "AXIN2"),
  gene_1_id = c(
    "ENSG00000133083", "ENSG00000164816", "ENSG00000122859",
    "ENSG00000100604", "ENSG00000198788", "ENSG00000139292"
  ),
  gene_2 = c(
    "ENSG00000095303", "ENSG00000164822", "ENSG00000173404",
    "ENSG00000157005", "ENSG00000275395", "ENSG00000168646"
  )
)

counts <- norm_counts %>%
  pivot_longer(starts_with("TCGA"), names_to="sample", values_to="exp") %>%
  pivot_wider(names_from="Ensembl_ID", values_from="exp") %>%
  transmute(
    sample=sample,
    DCLK1=ENSG00000133083, PTGS1=ENSG00000095303,
    DEFA5=ENSG00000164816, DEFA6=ENSG00000164822,
    NEUROG3=ENSG00000122859, INSM1=ENSG00000173404,
    SST=ENSG00000157005, CHGA=ENSG00000100604,
    FCGBP=ENSG00000275395, MUC2=ENSG00000198788,
    AXIN2=ENSG00000168646, LGR5=ENSG00000139292 
  )

log2_counts <- mutate(counts, across(-sample, ~log2(.x + 1)))

wes_colors <- as.character(wes_palette("Zissou1", 6, "continuous"))

i <- 1
plots <- pmap(comparisons, function(...) {
  args <- list(...)
  p <- ggplot(log2_counts, aes_string(x=args$gene_1_name, y=args$gene_2_name)) +
    geom_point(color=wes_colors[i]) +
    geom_smooth(method="lm", color="black", lty=2) +
    theme_classic() +
    theme(aspect.ratio=1)
  i <<- i + 1
  return(p)
})

plots <- wrap_plots(plots)

ggsave(
  file.path("results", "correlation_plot.pdf"), plot=plots,
  device=cairo_pdf, height=8, width=8
)

## Stats.

stats <- pmap_df(comparisons, function(...) {
  args <- list(...)
  results <- counts %>%
    mutate(across(-sample, ~log2(.x + 1))) %>%
    {cor.test(as.formula(str_c("~", args$gene_1_name, "+", args$gene_2_name)), data=.)}

  results <- data.frame(
    comparison=str_c(args$gene_1_name,"vs", args$gene_2_name, sep="_"),
    pval=results$p.value,
    pearson.r=results$estimate
  )
  return(results)
})

stats <- mutate(stats, padj=p.adjust(pval, method="BH"))

fwrite(
  stats, file.path("results", "correlation_stats.tsv"),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)
