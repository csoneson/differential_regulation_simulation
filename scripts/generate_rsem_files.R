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
print(nbr_diff_reg)
print(nbr_diff_spliced)
print(nbr_diff_expr)

## -----------------------------------------------------------------------------
## Read file with mean/dispersion relationship
## -----------------------------------------------------------------------------
meandisp <- readRDS(mean_disp_file)
meanvect <- meandisp$pickrell.cheung.mu
dispvect <- meandisp$pickrell.cheung.phi

## -----------------------------------------------------------------------------
## Help functions
## -----------------------------------------------------------------------------
rdirichlet2 <- function(n, alpha) {
    a <- rdirichlet(n, alpha)
    a[is.na(a)] <- 0
    a
}

## -----------------------------------------------------------------------------
## Define the isoform results file to start from (from RSEM). 
## We will modify the TPM column of this file to create the files that 
## are used in the simulation (the simulator uses only the TPM column, so 
## we will not change the others)
## -----------------------------------------------------------------------------
message("Reading isoform results file...")
isoform.initial <- read.delim(isoform_results_file, header = TRUE, as.is = TRUE)
    
message("Creating gene summary...")
gene_summary <- isoform.initial %>% group_by(gene_id) %>%
    summarise(expected_gene_count_gr1 = sum(expected_count),
              effective_gene_length = sum(effective_length * IsoPct/100),
              nbr_isoforms = length(IsoPct[!grepl("-U", transcript_id)]),
              nbr_unspliced_isoforms = length(IsoPct[grepl("-U", transcript_id)]), 
              nbr_expr_isoforms = length(which(IsoPct[!grepl("-U", transcript_id)] > 0)),
              nbr_expr_isoforms10 = length(which(IsoPct[!grepl("-U", transcript_id)] > 10)),
              sum_isopct_spliced_spliceable = sum(IsoPct[!grepl("-U", transcript_id) & 
                                                             (paste0(transcript_id, "-U") %in% 
                                                                  transcript_id)]),
              sum_isopct_unspliced = sum(IsoPct[grepl("-U", transcript_id)])) %>%
    dplyr::ungroup()

## -----------------------------------------------------------------------------
## Introduce differential expression
## -----------------------------------------------------------------------------
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr1
gene_summary$gene_de_status <- 0

## Generate fold changes if there are differentially expressed genes
if (nbr_diff_expr > 0) {
    set.seed(seed)
    fold_changes <- (2 + rexp(nbr_diff_expr, 
                              rate = 1))^(c(-1, 1)[round(runif(nbr_diff_expr)) + 1])
    
    diff_expr_genes <- sample(1:nrow(gene_summary), nbr_diff_expr, replace = FALSE)
    gene_summary$expected_gene_count_gr2[diff_expr_genes] <- 
        gene_summary$expected_gene_count_gr1[diff_expr_genes] * fold_changes
    gene_summary$gene_de_status[diff_expr_genes] <- 1
}
    
## -----------------------------------------------------------------------------
## Adjust the expected gene count to the desired library size, to obtain 
## the right dispersion estimates.
## Also calculate gene RPKs (proportional to TPMs)
## -----------------------------------------------------------------------------
gene_summary$expected_gene_count_gr1 <- gene_summary$expected_gene_count_gr1 /
    sum(gene_summary$expected_gene_count_gr1) * library_size
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr2 /
    sum(gene_summary$expected_gene_count_gr2) * library_size

gene_summary$dispersion_gr1 <- sapply(gene_summary$expected_gene_count_gr1, 
                                      function(w) {
                                          dispvect[which.min(abs(meanvect - w))]
                                      })
gene_summary$dispersion_gr2 <- sapply(gene_summary$expected_gene_count_gr2, 
                                      function(w) {
                                          dispvect[which.min(abs(meanvect - w))]
                                      })

message("Calculating gene counts and RPKs...")
for (i in 1:nbr_per_group) {
    gene_summary[, paste0("s", i, "_geneCount")] <- 
        sapply(1:nrow(gene_summary), function(j) {
            rnbinom(n = 1, size = 1/gene_summary$dispersion_gr1[j], 
                    mu = gene_summary$expected_gene_count_gr1[j])
        })
    gene_summary[, paste0("s", i, "_geneRPK")] <- 
        gene_summary[, paste0("s", i, "_geneCount")] /
        gene_summary$effective_gene_length * 1e3
}
for (i in (nbr_per_group + 1):(2 * nbr_per_group)) {
    gene_summary[, paste0("s", i, "_geneCount")] <- 
        sapply(1:nrow(gene_summary), function(j) {
            rnbinom(n = 1, size = 1/gene_summary$dispersion_gr2[j], 
                    mu = gene_summary$expected_gene_count_gr2[j])
        })
    gene_summary[, paste0("s", i, "_geneRPK")] <- 
        gene_summary[, paste0("s", i, "_geneCount")] /
        gene_summary$effective_gene_length * 1e3
}

## -----------------------------------------------------------------------------
## Create isoform summary
## -----------------------------------------------------------------------------
message("Creating isoform summary...")
isoform_summary <- isoform.initial
## Add some variability to the isoform percentages
set.seed(seed)
for (i in 1:(2 * nbr_per_group)) {
    isoform_summary <- isoform_summary %>% group_by(gene_id) %>% 
        mutate(IsoPctDirichlet = c(rdirichlet2(1, IsoPct/100 * 100))) %>%
        setNames(c(colnames(isoform_summary), paste0("s", i, "_IsoPct"))) %>%
        dplyr::ungroup()
}

## -----------------------------------------------------------------------------
## Introduce differential splicing (change the s*_IsoPct for the samples in group 2)
## -----------------------------------------------------------------------------
if (nbr_diff_spliced > 0) {
    ## Extract genes with at least 2 isoforms and expected count > 10
    message("Introducing differential splicing...")
    ds_genes <- 
        gene_summary$gene_id[sample(intersect(which(gene_summary$nbr_isoforms >= 2),
                                              which(gene_summary$expected_gene_count_gr1 > 10)), 
                                    nbr_diff_spliced, replace = FALSE)]
    isoform_summary_nonds <- isoform_summary[!(isoform_summary$gene_id %in% ds_genes), ]
    isoform_summary_ds <- isoform_summary[isoform_summary$gene_id %in% ds_genes, ]
    
    ## Helper function to modify isoform percentages
    mutate_IsoPct <- function(w, txid, getpctfrom) {
        ## w is a vector of isoform proportions
        ## txid is a vector of transcript IDs
        ## getpctfrom is a vector of transcript IDs to get new isoform pct from
        
        ## Find spliced transcript IDs
        spl <- txid[which(!grepl("-U", txid))]
        
        ## Get isoform percentages for spliced transcripts and corresponding 
        ## unspliced ones
        wspl <- w[match(spl, txid)]
        wuspl <- w[match(paste0(spl, "-U"), txid)]
        
        ## Calculate fraction of total isoform percentage that comes from the 
        ## spliced variant
        wtot <- rowSums(cbind(wspl, wuspl), na.rm = TRUE)
        relspl <- wspl/wtot
        
        ## Assign new total isoform percentage to transcripts
        wspl <- wtot[match(getpctfrom[match(spl, txid)], spl)]
        
        ## Split total isoform percentage between spliced and unspliced variant
        wuspl <- (1 - relspl) * wspl
        wspl <- relspl * wspl
        
        ## Return vector of new isoform percentages, in the same order as 
        ## the input
        names(wspl) <- spl
        names(wuspl) <- paste0(spl, "-U")
        wall <- c(wspl, wuspl)[txid]
        
        wall
    }
    
    ## Isoform percentages will be shuffled between transcripts - need to 
    ## make sure the shuffling is the same for all samples.
    get_shuffling <- function(txid) {
        splidx <- !grepl("-U", txid)
        tgt <- sample(txid[splidx], sum(splidx))
        while (all(tgt == txid[splidx])) {
            tgt <- sample(txid[splidx], sum(splidx))
        }
        map <- structure(rep(tgt, 2), names = c(txid[splidx], paste0(txid[splidx], "-U")))
        map[txid]
    }
    
    isoform_summary_ds <- isoform_summary_ds %>% group_by(gene_id) %>%
        mutate(getIsoPctFrom = get_shuffling(transcript_id)) %>%
        mutate(s4_IsoPct = mutate_IsoPct(s4_IsoPct, transcript_id, getIsoPctFrom),
               s5_IsoPct = mutate_IsoPct(s5_IsoPct, transcript_id, getIsoPctFrom),
               s6_IsoPct = mutate_IsoPct(s6_IsoPct, transcript_id, getIsoPctFrom),
               gene_ds_status = 1,
               transcript_ds_status = as.numeric(sub("-U", "", transcript_id) != getIsoPctFrom)) %>%
        dplyr::ungroup()
    
    isoform_summary_nonds <- isoform_summary_nonds %>% group_by(gene_id) %>%
        mutate(getIsoPctFrom = sub("-U", "", transcript_id)) %>%
        mutate(gene_ds_status = 0,
               transcript_ds_status = 0) %>%
        dplyr::ungroup()
    
    isoform_summary <- dplyr::bind_rows(isoform_summary_nonds, isoform_summary_ds)
} else {
    isoform_summary <- isoform_summary %>% group_by(gene_id) %>%
        mutate(getIsoPctFrom = sub("-U", "", transcript_id)) %>%
        mutate(gene_ds_status = 0,
               transcript_ds_status = 0) %>%
        dplyr::ungroup()
}

## -----------------------------------------------------------------------------
## Introduce differential regulation (swap the s*_IsoPct between -U and non-U 
## isoforms for the samples in group 2)
## -----------------------------------------------------------------------------
if (nbr_diff_reg > 0) {
    ## Extract genes with at least 1 unspliced isoform and expected count > 10
    message("Introducing differential regulation...")
    dr_genes <- 
        gene_summary$gene_id[sample(
            intersect(intersect(which(gene_summary$nbr_unspliced_isoforms >= 1),
                                which(gene_summary$expected_gene_count_gr1 > 10)),
                      which(abs(gene_summary$sum_isopct_spliced_spliceable / 
                                    (gene_summary$sum_isopct_spliced_spliceable + 
                                         gene_summary$sum_isopct_unspliced) - 0.5) > 0.1)
            ), 
            nbr_diff_reg, replace = FALSE)]
    isoform_summary_nondr <- isoform_summary[!(isoform_summary$gene_id %in% dr_genes), ]
    isoform_summary_dr <- isoform_summary[isoform_summary$gene_id %in% dr_genes, ]
    
    swap_IsoPct <- function(w, txid) {
        wtmp <- w
        for (i in grep("-U", txid)) {
            wtmp[i] <- w[txid == sub("-U", "", txid[i])]
            wtmp[txid == sub("-U", "", txid[i])] <- w[i]
        }
        wtmp
    }
    
    isoform_summary_dr <- isoform_summary_dr %>% group_by(gene_id) %>%
        mutate(s4_IsoPct = swap_IsoPct(s4_IsoPct, transcript_id),
               s5_IsoPct = swap_IsoPct(s5_IsoPct, transcript_id),
               s6_IsoPct = swap_IsoPct(s6_IsoPct, transcript_id),
               gene_dr_status = 1,
               transcript_dr_status = as.numeric(sub("-U-U", "-U", paste0(transcript_id, "-U")) %in%
                                                     transcript_id)) %>%
        dplyr::ungroup()
    
    isoform_summary_nondr <- isoform_summary_nondr %>% group_by(gene_id) %>%
        mutate(gene_dr_status = 0,
               transcript_dr_status = 0) %>%
        dplyr::ungroup()
    
    isoform_summary <- dplyr::bind_rows(isoform_summary_nondr, isoform_summary_dr)
} else {
    isoform_summary <- isoform_summary %>% group_by(gene_id) %>%
        mutate(gene_dr_status = 0,
               transcript_dr_status = 0) %>%
        dplyr::ungroup()
}

## -----------------------------------------------------------------------------
## Generate final summary table
## -----------------------------------------------------------------------------
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
        final_summary$effective_length / 1000
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
    final_summary[, i][!is.finite(final_summary[, i])] <- 0  ## happens when expected_count=0, TPM>0
    final_summary[, i] <- final_summary[, i]/sum(final_summary[, i]) * 1e6
}

## -----------------------------------------------------------------------------
## Write to files
## -----------------------------------------------------------------------------
message("Writing result files...")
write.table(final_summary, file = file.path(outdir_base, 
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
    write.table(tmp, file = file.path(outdir_base, "1_reads", 
                                      "rsem_files", paste0("sample", i, ".txt")), 
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

date()
sessionInfo()
