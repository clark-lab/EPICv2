mygrl <- function (peaks, TxDb, seqlev, nucleotideLevel = FALSE, ignore.strand = TRUE, 
         promoterRegion = c(upstream = 2000, downstream = 100), geneDownstream = c(upstream = 0, 
                                                                                   downstream = 1000), labels = list(geneLevel = c(promoter = "Promoter", 
                                                                                                                                   geneDownstream = "Downstream", geneBody = "Gene body", 
                                                                                                                                   distalIntergenic = "Distal Intergenic"), ExonIntron = c(exon = "Exon", 
                                                                                                                                                                                           intron = "Intron", intergenic = "Intergenic"), Exons = c(utr5 = "5' UTR", 
                                                                                                                                                                                                                                                    utr3 = "3' UTR", CDS = "CDS", otherExon = "Other exon"), 
                                                                                                                     group = c(geneLevel = "Gene Level", promoterLevel = "Promoter Level", 
                                                                                                                               Exons = "Exon level", ExonIntron = "Exon/Intron/Intergenic")), 
         labelColors = c(promoter = "#D55E00", geneDownstream = "#E69F00", 
                         geneBody = "#51C6E6", distalIntergenic = "#AAAAAA", exon = "#009DDA", 
                         intron = "#666666", intergenic = "#DDDDDD", utr5 = "#0072B2", 
                         utr3 = "#56B4E9", CDS = "#0033BF", otherExon = "#009E73"), 
         plot = TRUE, keepExonsInGenesOnly = TRUE, promoterLevel) 
{
  stopifnot(`peaks must be an object of GRanges or GRangesList` = inherits(peaks, 
                                                                           c("GRanges", "GRangesList")))
  if (is(peaks, "GRanges")) {
    n <- deparse(substitute(peaks))
    peaks <- GRangesList(peaks)
    names(peaks) <- n
    isGRanges <- TRUE
  }
  else {
    isGRanges <- FALSE
  }
  stopifnot(`TxDb must be an object of TxDb` = is(TxDb, "TxDb"))
  stopifnot(`nuleotideLevel is not logical` = is.logical(nucleotideLevel))
  stopifnot(`promoterRegion should contain element upstream and downstream` = all(c("upstream", 
                                                                                    "downstream") %in% names(promoterRegion)))
  stopifnot(`geneDownstream should contain element upstream and downstream` = all(c("upstream", 
                                                                                    "downstream") %in% names(geneDownstream)))
  stopifnot(`Elements in promoterRegion should be numeric` = is.numeric(promoterRegion))
  stopifnot(`Elements in geneDownstream should be numeric` = is.numeric(geneDownstream))
  labs <- list(geneLevel = c(promoter = "Promoter", geneDownstream = "Downstream", 
                             geneBody = "Gene body", distalIntergenic = "Distal Intergenic"), 
               ExonIntron = c(exon = "Exon", intron = "Intron", intergenic = "Intergenic"), 
               Exons = c(utr5 = "5' UTR", utr3 = "3' UTR", CDS = "CDS", 
                         otherExon = "Other exon"))
  groupLabels <- c(geneLevel = "Gene Level", promoterLevel = "Promoter Level", 
                   Exons = "Exon level", ExonIntron = "Exon/Intron/Intergenic")
  labelCols = c(promoter = "#D55E00", geneDownstream = "#E69F00", 
                geneBody = "#51C6E6", distalIntergenic = "#AAAAAA", exon = "#009DDA", 
                intron = "#666666", intergenic = "#DDDDDD", utr5 = "#0072B2", 
                utr3 = "#56B4E9", CDS = "#0033BF", otherExon = "#009E73", 
                undefined = "#FFFFFF")
  labelCols[names(labelColors)] <- labelColors
  for (i in names(labs)) {
    if (i %in% names(labels)) {
      n <- intersect(names(labs[[i]]), names(labels[[i]]))
      labs[[i]] <- c(labels[[i]][n], labs[[i]][!names(labs[[i]]) %in% 
                                                 n])
    }
  }
  for (i in names(groupLabels)) {
    if ("group" %in% names(labels)) {
      groupLabels[names(labels[["group"]])] <- labels[["group"]]
    }
  }
  if (!missing(promoterLevel)) {
    promoterLevel$breaks <- sort(promoterLevel$breaks)
    stopifnot(`breaks, labels and colors of promoterLevel are not paired` = length(promoterLevel$breaks) == 
                length(promoterLevel$labels) + 1 && length(promoterLevel$labels) == 
                length(promoterLevel$colors))
    proK <- paste0("promoter", seq_along(promoterLevel$labels))
    promoterLevel$upstream <- ifelse(promoterLevel$breaks < 
                                       0, abs(promoterLevel$breaks), 0)
    promoterLevel$upstream <- promoterLevel$upstream[-length(promoterLevel$upstream)]
    promoterLevel$downstream <- ifelse(promoterLevel$breaks > 
                                         0, promoterLevel$breaks, 0)
    promoterLevel$downstream <- promoterLevel$downstream[-1]
    proV <- promoterLevel$labels
    names(proV) <- proK
    labs <- c(list("promoterLevel" = proV), labs)
    proV <- promoterLevel$colors
    names(proV) <- proK
    labelCols <- c(proV, labelCols)
  }
  else {
    promoterLevel <- NULL
  }
  defaultW <- getOption("warn")
  options(warn = -1)
  on.exit({
    warn = defaultW
  })
  seql <- seqlevelsStyle(peaks)
  anno <- GRangesList()
  for (i in names(labs)) {
    anno[[i]] <- switch(i, "promoterLevel" = {
      suppressMessages(g <- genes(TxDb, single.strand.genes.only = TRUE))
      pro <- mapply(FUN = function(upstream, downstream) {
        promoters(g, upstream = upstream, downstream = downstream)
      }, promoterLevel$upstream, promoterLevel$downstream, 
      SIMPLIFY = FALSE)
      names(pro) <- names(labs[["promoterLevel"]])
      pro <- rev(pro)
      current_anno <- GRanges()
      for (j in seq_along(pro)) {
        ca <- ChIPpeakAnno:::filterByOverlaps(pro[[j]], current_anno, 
                               ignore.strand = ignore.strand)
        mcols(ca) <- DataFrame(type = rep(names(pro)[j], 
                                          length(ca)))
        current_anno <- c(current_anno, ca)
      }
      seqlevelsStyle(current_anno) <- seql[1]
      current_anno
    }, "geneLevel" = {
      suppressMessages(g <- genes(TxDb, single.strand.genes.only = TRUE))
      pro <- promoters(g, upstream = promoterRegion["upstream"], 
                       downstream = promoterRegion["downstream"])
      dws <- downstreams(g, upstream = geneDownstream["upstream"], 
                         downstream = geneDownstream["downstream"])
      pro <- GenomicRanges::trim(pro)
      dws <- GenomicRanges::trim(dws)
      intergenic <- gaps(GenomicRanges::reduce(c(pro, g, dws), ignore.strand = FALSE))
      intergenic <- intergenic[!strand(intergenic) %in% "*"]
      current_anno <- GRanges()
      for (j in names(labs[["geneLevel"]])) {
        s <- switch(j, "promoter" = pro, "geneDownstream" = dws, 
                    "geneBody" = g, "distalIntergenic" = intergenic)
        ca <- ChIPpeakAnno:::filterByOverlaps(s, current_anno, ignore.strand = ignore.strand)
        mcols(ca) <- DataFrame(type = rep(j, length(ca)))
        current_anno <- c(current_anno, ca)
      }
      seqlevelsStyle(current_anno) <- seql[1]
      current_anno
    }, "ExonIntron" = {
      exon <- exons(TxDb)
      intron <- unlist(intronsByTranscript(TxDb))
      if (keepExonsInGenesOnly) {
        suppressMessages(g <- genes(TxDb, single.strand.genes.only = TRUE))
        ole <- findOverlaps(exon, g, type = "within")
        oli <- findOverlaps(intron, g, type = "within")
        ole <- !seq_along(exon) %in% queryHits(ole)
        oli <- !seq_along(intron) %in% queryHits(oli)
        if (sum(ole) > 0 || sum(oli) > 0) {
          warning(paste(sum(ole), "exons were dropped because there is no", 
                        "relative gene level annotations.", sum(oli), 
                        "introns were dropped because there is no", 
                        "relative gene level annotations."))
          exon <- exon[!ole]
          intron <- intron[!oli]
        }
      }
      intergenic <- gaps(GenomicRanges::reduce(c(exon, intron), ignore.strand = FALSE))
      intergenic <- intergenic[!strand(intergenic) %in% "*"]
      current_anno <- GRanges()
      for (j in names(labs[["ExonIntron"]])) {
        s <- switch(j, "exon" = exon, "intron" = intron, 
                    "intergenic" = intergenic)
        ca <- ChIPpeakAnno:::filterByOverlaps(s, current_anno, ignore.strand = ignore.strand)
        mcols(ca) <- DataFrame(type = rep(j, length(ca)))
        current_anno <- c(current_anno, ca)
      }
      seqlevelsStyle(current_anno) <- seql[1]
      current_anno
    }, "Exons" = {
      utr5 <- unlist(fiveUTRsByTranscript(TxDb))
      utr3 <- unlist(threeUTRsByTranscript(TxDb))
      CDS <- cds(TxDb)
      exon <- exons(TxDb)
      if (keepExonsInGenesOnly) {
        suppressMessages(g <- genes(TxDb, single.strand.genes.only = TRUE))
        ole <- findOverlaps(exon, g, type = "within")
        olc <- findOverlaps(CDS, g, type = "within")
        ol5 <- findOverlaps(utr5, g, type = "within")
        ol3 <- findOverlaps(utr3, g, type = "within")
        ole <- !seq_along(exon) %in% queryHits(ole)
        olc <- !seq_along(CDS) %in% queryHits(olc)
        ol5 <- !seq_along(utr5) %in% queryHits(ol5)
        ol3 <- !seq_along(utr3) %in% queryHits(ol3)
        if (sum(ole) > 0 || sum(olc) > 0 || sum(ol5) > 
            0 || sum(ol3)) {
          warning(paste(sum(ole), "exons were dropped because there is no", 
                        "relative gene level annotations.", sum(olc), 
                        "CDS were dropped because there is no", "relative gene level annotations.", 
                        sum(ol5), "utr5 were dropped because there is no", 
                        "relative gene level annotations.", sum(ol3), 
                        "utr3 were dropped because there is no", 
                        "relative gene level annotations."))
          exon <- exon[!ole]
          CDS <- CDS[!olc]
          utr5 <- utr5[!ol5]
          utr3 <- utr3[!ol3]
        }
      }
      current_anno <- GRanges()
      for (j in names(labs[["Exons"]])) {
        s <- switch(j, "utr5" = utr5, "utr3" = utr3, "CDS" = CDS, 
                    "otherExon" = exon)
        ca <- ChIPpeakAnno:::filterByOverlaps(s, current_anno, ignore.strand = ignore.strand)
        mcols(ca) <- DataFrame(type = rep(j, length(ca)))
        current_anno <- c(current_anno, ca)
      }
      seqlevelsStyle(current_anno) <- seql[1]
      current_anno
    })
  }
  groupLabels <- groupLabels[names(anno)]
  groupLabels[is.na(groupLabels)] <- names(anno)[is.na(groupLabels)]
  if (!missing(seqlev)) {
    if (length(seqlev) > 0) {
      peaks <- lapply(peaks, function(.ele) .ele[seqnames(.ele) %in% 
                                                   seqlev])
    }
  }
  if (nucleotideLevel) {
    peaks <- lapply(peaks, function(.ele) {
      y <- disjoin(c(.ele, unlist(anno)), ignore.strand = ignore.strand)
      subsetByOverlaps(y, .ele, ignore.strand = ignore.strand)
    })
  }
  peaks <- lapply(peaks, FUN = function(.peaks) {
    pct <- lapply(anno, FUN = function(.ele) {
      y <- .peaks
      ol <- findOverlaps(y, .ele, ignore.strand = ignore.strand)
      ol <- as.data.frame(ol)
      ol <- ol[order(ol$queryHits, ol$subjectHits), ]
      ol <- ol[!duplicated((ol$queryHits)), ]
      y$anno <- rep("undefined", length(y))
      y$anno[ol$queryHits] <- .ele$type[ol$subjectHits]
      y$anno
    })
    pct <- do.call(cbind, pct)
    mcols(.peaks) <- cbind(mcols(.peaks), pct)
    .peaks
  })
  if (isGRanges) {
    peaks <- peaks[[1]]
  }
  melt <- function(m) {
    if (nucleotideLevel) {
      m <- m[rep(seq_along(m), width(m))]
    }
    m <- mcols(m)[, names(anno)]
    p <- lapply(names(anno), function(.ele) {
      tt <- table(m[, .ele])
      data.frame(category = factor(rep(groupLabels[.ele], 
                                       length(tt)), levels = groupLabels), type = factor(names(tt), 
                                                                                         levels = rev(names(labelCols))), percentage = as.numeric(tt)/sum(tt))
    })
    do.call(rbind, p)
  }
  if (isGRanges) {
    dat <- melt(peaks)
    l <- unlist(unname(labs))
    l1 <- paste0(l, " (", round(dat$percentage[match(names(l), 
                                                     dat$type)] * 100, digits = 1), "%)")
    names(l1) <- names(l)
    l1 <- c(l1, undefined = "")
    p <- ggplot(dat, aes_string(x = "category", y = "percentage", 
                                fill = "type")) + geom_col() + coord_polar("y") + 
      geom_text(data = subset(dat, !duplicated(dat$category)), 
                aes_string(x = "category", label = "category"), 
                y = 1) + theme_void()
  }
  else {
    dat <- lapply(peaks, melt)
    dat1 <- do.call(rbind, dat)
    dat1$source <- rep(names(peaks), vapply(dat, nrow, FUN.VALUE = 0))
    l1 <- c(unlist(unname(labs)), undefined = "")
    p <- ggplot(dat1, aes_string(x = "source", y = "percentage", 
                                 fill = "type")) + geom_bar(stat = "identity", position=position_dodge()) + coord_flip() + 
      facet_wrap(as.formula("~ category"), ncol = 1) + 
      theme_bw()
  }
  p <- p + scale_fill_manual(values = labelCols, labels = l1, 
                             name = NULL, guide = guide_legend(reverse = TRUE))
  if (plot) {
    print(p)
  }
  return(invisible(list(peaks = peaks, plot = p)))
  return(dat1)
}