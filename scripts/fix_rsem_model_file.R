args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(inputFile)
print(outputFile)

inp <- readLines(inputFile)

## The first 13 lines are not modified (model type, forward probability, 
## fragment length and read start position distribution)
## Next, the quality score distribution is modified
## Line 14 is the initial probabilities of each quality score
## Lines 15-114 are the transition probabilities
for (i in 14:114) {
    tmp <- as.numeric(strsplit(inp[i], " ")[[1]])
    if (any(tmp != 0)) {
        ## Set value for quality score 2 to 0 and rescale the rest to sum to 1
        tmp[3] <- 0
        tmp <- tmp/sum(tmp)
        inp[i] <- paste(tmp, collapse = " ")
    }
}

## Next is the sequencing error model, in 100 blocks of 5 rows.
## Set all values in the third block (quality = 2) to 0
## line 115 empty, line 116 summary info, block 3 from line 129-133
for (i in 129:133) {
    inp[i] <- "0 0 0 0 0"
}
## Last block ends at line 116 + (100 * 6) - 1 = 715
## Next there's an empty line, a summary line, and again we want to set the 
## third line after that (quality score 2) to all zeros
inp[720] <- "0 0 0 0 0"

## Save output file
writeLines(inp, con = outputFile)

date()
sessionInfo()
