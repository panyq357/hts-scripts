library(ggplot2)
library(patchwork)
library(GenomicRanges)

main <- function() {

    gfp_bw <- rtracklayer::import("chipseq-pipeline/results/bamCoverage/GFP.bw")
    pc48_bw <- rtracklayer::import("chipseq-pipeline/results/bamCoverage/PC48.bw")
    pc54_bw <- rtracklayer::import("chipseq-pipeline/results/bamCoverage/PC54.bw")

    pc48_peak <- rtracklayer::import("chipseq-pipeline/results/macs_callpeak/PC48/PC48_peaks.narrowPeak")
    pc54_peak <- rtracklayer::import("chipseq-pipeline/results/macs_callpeak/PC54/PC54_peaks.narrowPeak")

    gtf <- rtracklayer::import("/home/panyq/Tools/index-scripts/os/rap-db/results/gtf/os.rap-db.make.gtf")

    genes <- readxl::read_excel("Genes.xlsx")

    plot_all <- function(gene, name) {

        gene_range <- gene_to_range(gtf, gene, upstream=2000, downstream=2000)

        svg(sprintf("results/plot_track/%s.svg", name), width=5, height=5)
        print(
            (plot_bw(gfp_bw, gene_range, ylim=c(0, 15), name="GFP") + theme(axis.text.x=element_blank())) +
            (plot_bw(pc48_bw, gene_range, ylim=c(0, 15), name="PC48") + theme(axis.text.x=element_blank())) +
            (plot_bed(pc48_peak, gene_range, name="PC48\nPeak") + theme(axis.text.x=element_blank())) +
            (plot_bw(pc54_bw, gene_range, ylim=c(0, 15), name="PC54") + theme(axis.text.x=element_blank())) +
            (plot_bed(pc54_peak, gene_range, name="PC54\nPeak") + theme(axis.text.x=element_blank())) +
            (plot_gene_track(gtf, gene_range, keep_only=gene) + labs(y=name) + theme(axis.title.y=element_text(angle=0, vjust=0.5))) +
            plot_layout(design = "1\n2\n3\n4\n5\n6", heights = c(3, 3, 0.5, 3, 0.5, 3))
        )
        dev.off()
    }

    if (!dir.exists("results/plot_track")) dir.create("results/plot_track")

    for (i in 1:nrow(genes)) {
        plot_all(genes[[i, "ID"]], genes[[i, "Name"]])
    }
}

#' Given a gene name, return its nearby region GRanges.
gene_to_range <- function(gtf, gene, upstream=2000, downstream=2000) {
    x <- subset(gtf, gene_id == gene) |> range()
    start(x) <- start(x) - upstream
    end(x) <- end(x) + downstream
    strand(x) <- "*"
    return(x)
}

#' Given a GRanges of one gene, find its primary transcript, and return its GRanges.
get_primary_tx <- function(gene_gtf, by="cds") {
    tx_list <- split(gene_gtf, gene_gtf$transcript_id)

    if (by == "cds") {
        x <- sapply(tx_list, function(tx) {
            subset(tx, type == "CDS") |> width() |> sum()
        }) |> sort(decreasing=T)
        primary_tx_id <- names(x)[[1]]
    }

    primary_tx <- tx_list[[primary_tx_id]]
    return(primary_tx)
}

#' Given a list of GRanges, distribute their ranges into different tracks
get_track_list <- function(grl) {
    range_list <- lapply(grl, range) |> lapply(unstrand)
    track_list <- list(setNames(list(range_list[[1]]), names(range_list)[[1]]))

    if (length(range_list) == 1) {
        return(track_list)
    }

    for (i in 2:length(range_list)) {

        this_range <- range_list[[i]]
        this_range_name <- names(range_list)[[i]]

        j <- 0
        overlap <- TRUE
        while (overlap == TRUE) {
            j <- j + 1
            if (j > length(track_list)) {
                track_list[[j]] <- list()
                overlap <- FALSE
            } else {
                overlap <- any(sapply(track_list[[j]], function(x) x %over% this_range))
            }
        }
        track_list[[j]][[this_range_name]] <- this_range
    }

    return(track_list)
}

add_track_num <- function(grl) {
    track_list <- get_track_list(grl)

    track_num_list <- rev(seq(2, 2 * length(track_list), length.out=length(track_list)))

    for (i in 1:length(track_list)) {
        for (name in names(track_list[[i]])) {
            grl[[name]]$track_num <- track_num_list[[i]]
        }
    }
    return(grl)
}

get_arrow_pos <- function(plot_data, width=500) {
    cds <- subset(plot_data, type == "CDS")
    if (all(cds$strand == "+")) {
        start <- min(cds$start)
        end <- start + width
    } else if (all(cds$strand == "-")) {
        start <- max(cds$end)
        end <- start - width
    }
    return(c(start=start, end=end, track_num=cds$track_num[[1]]))
}


plot_gene_track <- function(gtf, gene_range, remove=FALSE, keep_only=FALSE) {

    gene_in_range <- unique(subsetByOverlaps(gtf, gene_range)$gene_id)

    if (remove !=FALSE) {
        gene_in_range <- gene_in_range[!gene_in_range %in% remove]
    }

    if (keep_only != FALSE) {
        gene_in_range <- keep_only
    }

    range_gr <- subset(gtf, gene_id %in% gene_in_range)

    plot_data <- split(range_gr, ~gene_id) |> 
        lapply(get_primary_tx) |>
        add_track_num() |>
        lapply(as.data.frame) |>
        do.call(what=rbind, args=_)

    arrow_data <- subset(plot_data, type=="CDS") |>
        split(~transcript_id) |>
        lapply(get_arrow_pos) |>
        do.call(what=rbind, args=_) |>
        as.data.frame()


    track_plot <- ggplot() +
        geom_rect(data=subset(plot_data, type=="transcript"), mapping=aes(xmin=start, xmax=end, ymin=track_num-0.05, ymax=track_num+0.05), fill="black", color="black") +
        geom_rect(data=subset(plot_data, type=="exon"), mapping=aes(xmin=start, xmax=end, ymin=track_num-0.5, ymax=track_num+0.5), fill="white", color="black") +
        geom_rect(data=subset(plot_data, type=="CDS"), mapping=aes(xmin=start, xmax=end, ymin=track_num-0.5, ymax=track_num+0.5), fill="black", color="black") +
        geom_text(data=subset(plot_data, type=="transcript"), mapping=aes(x=(start+end)/2, y=track_num-1, label=gene_id)) +
        geom_segment(data=arrow_data, mapping=aes(x=start, y=track_num+1, xend=end, yend=track_num+1), arrow=arrow(angle=45, length=unit(0.1, "inches"))) +
        geom_segment(data=arrow_data, mapping=aes(x=start, y=track_num, xend=start, yend=track_num+1)) +
        coord_cartesian(xlim=c(start(gene_range), end(gene_range)), ylim=c(min(plot_data$track_num)-2, max(plot_data$track_num)+2)) +
        scale_x_continuous(expand=c(0, 1)) +
        scale_y_continuous(breaks=NULL, labels=NULL) +
        labs(x=NULL, y=NULL) +
        theme_classic()
    
    return(track_plot)
}

plot_bw <- function(bw, gene_range, ylim, name="") {

    bw_plot <- subsetByOverlaps(bw, gene_range) |>
        as.data.frame() |>
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = score), fill="grey", color="grey", linewidth=0.1) +
        coord_cartesian(xlim=c(start(gene_range), end(gene_range)), ylim=ylim) +
        scale_x_continuous(expand=c(0, 1)) +
        labs(y=name) +
        theme_classic()

    return(bw_plot)
}

plot_bed <- function(bed, gene_range, name="") {
    
    bed_plot <- subsetByOverlaps(bed, gene_range) |>
        as.data.frame() |>
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0.5, ymax = 1.5), fill="black") +
        coord_cartesian(xlim=c(start(gene_range), end(gene_range)), ylim=c(0, 2)) +
        scale_x_continuous(expand=c(0, 1)) +
        scale_y_continuous(breaks=NULL, labels=NULL) +
        labs(y=name) +
        theme_classic()

    return(bed_plot)
}


main()
