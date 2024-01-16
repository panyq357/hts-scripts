#' Export normalized counts to a GCT file.
#' GCT format ref: <https://docs.gsea-msigdb.org/#GSEA/Data_Formats/#gct-gene-cluster-text-file-format-gct>
#' nc: normalized counts
write_gct <- function(nc, out) {
    lines = c(
        "#1.2",
        sprintf("%s\t%s", nrow(nc), ncol(nc)),
        paste(c("NAME", "Description", colnames(nc)), collapse="\t")
    )
    writeLines(lines, out)
    write.table(cbind(rownames(nc), NA, nc), out, row.names=F, col.names=F, append=T, quote=F, sep="\t")
}

#' Export coldata to a CLS file.
#' CLS format ref: <https://docs.gsea-msigdb.org/#GSEA/Data_Formats/#cls-categorical-eg-tumor-vs-normal-class-file-format-cls>
write_cls <- function(coldata, out) {
    lines = c(
        sprintf("%s %s 1", nrow(coldata), length(unique(coldata$Group))),
        sprintf("# %s", paste(unique(coldata$Group), collapse=" ")),
        sprintf(paste(coldata$Group, collapse=" "))
    )
    writeLines(lines, out)
}

