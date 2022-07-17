args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(gtools)
    library(dplyr)
})

print(seed)
print(isoform_results_file)
print(nbr_per_group)
print(mean_disp_file)
print(outdir_base)
print(library_size)
print(nbr_diff_spliced)
print(nbr_diff_expr)

## Generate fold changes if there are differentially expressed genes
if (nbr_diff_expr > 0) {
    set.seed(seed)
    fold_changes <- (2 + rexp(nbr_diff_expr, 
                             rate = 1))^(c(-1, 1)[round(runif(nbr_diff_expr)) + 1])
}

print(head(fold_changes))

## Read file with mean/dispersion relationship
meandisp <- readRDS(mean_disp_file)
meanvect <- meandisp$pickrell.cheung.mu
dispvect <- meandisp$pickrell.cheung.phi

## Generate RSEM files for individual samples, starting from the isoform 
## file estimated by RSEM from the real data
set.seed(seed)
    
rdirichlet2 <- function(n, alpha) {
    a <- rdirichlet(n, alpha)
    a[is.na(a)] <- 0
    a
}

## Define the isoform results file to start from (from RSEM). 
## We will modify the TPM column of this file to create the files that 
## are used in the simulation (the simulator uses only the TPM column, so 
## we will not change the others)
message("Reading isoform results file...")
isoform.initial <- read.delim(isoform_results_file, header = TRUE, as.is = TRUE)
    
message("Creating gene summary...")
gene_summary <- isoform.initial %>% group_by(gene_id) %>%
    summarise(expected_gene_count_gr1 = sum(expected_count),
              effective_gene_length = sum(effective_length * IsoPct/100),
              nbr_isoforms = length(IsoPct),
              nbr_expr_isoforms = length(which(IsoPct > 0)),
              nbr_expr_isoforms10 = length(which(IsoPct > 10)))

## Introduce differential expression
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr1
gene_summary$gene_de_status <- 0
if (nbr_diff_expr > 0) {
    diff_expr_genes <- sample(1:nrow(gene_summary), nbr_diff_expr, replace = FALSE)
    gene_summary$expected_gene_count_gr2[diff_expr_genes] <- 
        gene_summary$expected_gene_count_gr1[diff_expr_genes] * fold_changes
    gene_summary$gene_de_status[diff_expr_genes] <- 1
}
    
## Adjust the expected gene count to the desired library size, to obtain 
## the right dispersion estimates
gene_summary$expected_gene_count_gr1 <- gene_summary$expected_gene_count_gr1/
    sum(gene_summary$expected_gene_count_gr1) * library_size
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr2/
    sum(gene_summary$expected_gene_count_gr2) * library_size

gene_summary$dispersion_gr1 <- sapply(gene_summary$expected_gene_count_gr1, 
                                      function(w) {
                                          dispvect[which.min(abs(meanvect - w))]
                                      })
gene_summary$dispersion_gr2 <- sapply(gene_summary$expected_gene_count_gr2, 
                                      function(w) {
                                          dispvect[which.min(abs(meanvect - w))]
                                      })

message("Calculating gene counts...")
for (i in 1:nbr_per_group) {
    gene_summary[, paste0("s", i, "_geneCount")] <- 
        sapply(1:nrow(gene_summary), function(j) {
            rnbinom(n = 1, size = 1/gene_summary$dispersion_gr1[j], 
                    mu = gene_summary$expected_gene_count_gr1[j])
        })
    gene_summary[, paste0("s", i, "_geneRPK")] <- 
        gene_summary[, paste0("s", i, "_geneCount")]/
        gene_summary$effective_gene_length * 1e3
}
for (i in (nbr_per_group + 1):(2 * nbr_per_group)) {
    gene_summary[, paste0("s", i, "_geneCount")] <- 
        sapply(1:nrow(gene_summary), function(j) {
            rnbinom(n = 1, size = 1/gene_summary$dispersion_gr2[j], 
                    mu = gene_summary$expected_gene_count_gr2[j])
        })
    gene_summary[, paste0("s", i, "_geneRPK")] <- 
        gene_summary[, paste0("s", i, "_geneCount")]/
        gene_summary$effective_gene_length * 1e3
}

message("Creating isoform summary...")
isoform_summary <- isoform.initial
## Add some variability to the isoform percentages
for (i in 1:(2 * nbr_per_group)) {
    isoform_summary <- isoform_summary %>% group_by(gene_id) %>% 
        mutate(IsoPctDirichlet = c(rdirichlet2(1, IsoPct/100 * 100))) %>%
        setNames(c(colnames(isoform_summary), paste0("s", i, "_IsoPct"))) %>%
        ungroup()
}

isoform_summary_tmp <- isoform_summary %>% group_by(gene_id) %>% 
    mutate(diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1]/100,
           gene_ds_status = 0,
           transcript_ds_status = 0)

## Null simulation, final summary table
message("Merging gene and isoform summaries...")
final_summary <- merge(isoform_summary_tmp, gene_summary, by.x = "gene_id", 
                       by.y = "gene_id", all = TRUE)
for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformRPK")] <- 
        final_summary[, paste0("s", i, "_geneRPK")] * 
        final_summary[, paste0("s", i, "_IsoPct")]
}

for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformCount")] <- 
        final_summary[, paste0("s", i, "_isoformRPK")] * 
        final_summary$effective_length/1000
}

for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformTPM")] <- 
        round(final_summary[, paste0("s", i, "_isoformCount")] / 
                  final_summary$expected_count * final_summary$TPM, 2)
}

final_summary[is.na(final_summary)] <- 0

## Scale each isoformTPM column so that it sums to 1 million
idx <- grep("_isoformTPM", colnames(final_summary))
for (i in idx) {
    final_summary[, i] <- final_summary[, i]/sum(final_summary[, i]) * 1e6
}


## Write to files
message("Writing result files...")
write.table(final_summary, file = file.path(outdir_base, "null_simulation", "3_truth", 
                                            "simulation_details.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

for (i in 1:(2 * nbr_per_group)) {
    tmp <- isoform.initial
    tmp$TPM <- final_summary[match(tmp$transcript_id, 
                                   final_summary$transcript_id), 
                             paste0("s", i, "_isoformTPM")]
    tmp$TPM <- round(tmp$TPM, digits = 2)
    tmp$TPM <- as.character(tmp$TPM)
    tmp$TPM[tmp$TPM == "0"] <- "0.00"
    tmp$expected_count <- as.character(tmp$expected_count)
    tmp$expected_count[tmp$expected_count == "0"] <- "0.00"
    tmp$FPKM <- as.character(tmp$FPKM)
    tmp$FPKM[tmp$FPKM == "0"] <- "0.00"
    tmp$IsoPct <- as.character(tmp$IsoPct)
    tmp$IsoPct[tmp$IsoPct == "0"] <- "0.00"
    write.table(tmp, file = file.path(outdir_base, "null_simulation", "1_reads", 
                                      "rsem_files", paste0("sample", i, ".txt")), 
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}


## Introduce differential splicing (change the s*_IsoPct for the samples in group 2)
## and generate new final_summary
## Extract 1,000 genes with at least 2 isoforms and at least 500 expected counts per million
## in the RSEM data, and at least two isoforms with > 10% relative abundance
message("Introducing differential splicing...")
ds_genes <- 
    gene_summary$gene_id[sample(intersect(which(gene_summary$nbr_expr_isoforms10 >= 2),
                                          which(gene_summary$expected_gene_count_gr1 > 500)), 
                                nbr_diff_spliced, replace = FALSE)]
isoform_summary_nonds <- isoform_summary[!(isoform_summary$gene_id %in% ds_genes), ]
isoform_summary_ds <- isoform_summary[isoform_summary$gene_id %in% ds_genes, ]

mutate_IsoPct <- function(w, ref) {
    o <- order(ref, decreasing = TRUE)[1:2]
    w[rev(o)] <- w[o]
    w
}
mutated_IsoPct <- function(ref) {
    o <- order(ref, decreasing = TRUE)[1:2]
    w <- rep(0, length(ref))
    w[o] <- 1
    w
}

isoform_summary_ds <- isoform_summary_ds %>% group_by(gene_id) %>%
    mutate(s4_IsoPct = mutate_IsoPct(s4_IsoPct, IsoPct),
           s5_IsoPct = mutate_IsoPct(s5_IsoPct, IsoPct),
           s6_IsoPct = mutate_IsoPct(s6_IsoPct, IsoPct),
           diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1]/100,
           gene_ds_status = 1,
           transcript_ds_status = mutated_IsoPct(IsoPct))

isoform_summary_nonds <- isoform_summary_nonds %>% group_by(gene_id) %>%
    mutate(diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1]/100,
           gene_ds_status = 0,
           transcript_ds_status = 0)

isoform_summary <- rbind(isoform_summary_nonds, isoform_summary_ds)

## Non-null simulation, final summary table
message("Merging gene and isoform summaries...")
final_summary <- merge(isoform_summary, gene_summary, by.x = "gene_id", 
                       by.y = "gene_id", all = TRUE)
for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformRPK")] <- 
        final_summary[, paste0("s", i, "_geneRPK")] * 
        final_summary[, paste0("s", i, "_IsoPct")]
}

for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformCount")] <- 
        final_summary[, paste0("s", i, "_isoformRPK")] * 
        final_summary$effective_length/1000
}

for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformTPM")] <- 
        round(final_summary[, paste0("s", i, "_isoformCount")] / 
                  final_summary$expected_count * final_summary$TPM, 2)
}

final_summary[is.na(final_summary)] <- 0

## Scale each isoformTPM column so that it sums to 1 million
idx <- grep("isoformTPM", colnames(final_summary))
for (i in idx) {
    final_summary[, i] <- final_summary[, i]/sum(final_summary[, i]) * 1e6
}

## Write to files
message("Writing result files...")
write.table(final_summary, file = file.path(outdir_base, "non_null_simulation", 
                                            "3_truth", "simulation_details.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

for (i in 1:(2 * nbr_per_group)) {
    tmp <- isoform.initial
    tmp$TPM <- final_summary[match(tmp$transcript_id, 
                                   final_summary$transcript_id), 
                             paste0("s", i, "_isoformTPM")]
    tmp$TPM <- round(tmp$TPM, digits = 2)
    tmp$TPM <- as.character(tmp$TPM)
    tmp$TPM[tmp$TPM == "0"] <- "0.00"
    tmp$expected_count <- as.character(tmp$expected_count)
    tmp$expected_count[tmp$expected_count == "0"] <- "0.00"
    tmp$FPKM <- as.character(tmp$FPKM)
    tmp$FPKM[tmp$FPKM == "0"] <- "0.00"
    tmp$IsoPct <- as.character(tmp$IsoPct)
    tmp$IsoPct[tmp$IsoPct == "0"] <- "0.00"
    write.table(tmp, file = file.path(outdir_base, "non_null_simulation", "1_reads", 
                                      "rsem_files", paste0("sample", i, ".txt")), 
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

date()
sessionInfo()
