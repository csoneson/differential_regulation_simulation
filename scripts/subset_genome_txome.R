args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(genomein)
print(genomeout)
print(txomein)
print(txomeout)
print(chrnameout)
print(gtfin)
tmp <- strsplit(coordkeep, "_")[[1]]
chrkeep <- tmp[1]
posstart <- as.numeric(tmp[2])
posend <- as.numeric(tmp[3])
print(chrkeep)
print(posstart)
print(posend)

suppressPackageStartupMessages({
    library(Biostrings)
    library(rtracklayer)
})

genomein <- readDNAStringSet(genomein)
txomein <- readDNAStringSet(txomein)
gtfin <- rtracklayer::import(gtfin)

genome <- genomein[sapply(strsplit(names(genomein), " "), .subset, 1) %in% chrkeep]
genes <- gtfin$gene_id[as.vector(gtfin$type == "gene" & 
                                     seqnames(gtfin) == chrkeep & 
                                     start(gtfin) > posstart & 
                                     end(gtfin) < posend)]
txome <- txomein[names(txomein) %in% 
                     gtfin$transcript_id[as.vector(gtfin$type == "transcript" & 
                                                       gtfin$gene_id %in% genes)], ]
print(table(seqnames(gtfin[gtfin$transcript_id %in% names(txome)])))

writeXStringSet(genome, filepath = genomeout)
writeXStringSet(txome, filepath = txomeout)
write.table(chrkeep, file = chrnameout, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

date()
sessionInfo()
