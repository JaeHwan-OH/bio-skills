#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

get_arg <- function(key, default=NULL) {
  k <- which(args == key)
  if (length(k) == 0) return(default)
  if (k == length(args)) return(default)
  return(args[k+1])
}

rds <- get_arg("--rds")
out <- get_arg("--out")
genome <- get_arg("--genome", "mm10")
export_rna <- tolower(get_arg("--export_rna", "no"))
export_atac <- tolower(get_arg("--export_atac", "yes"))
fragments_override <- get_arg("--fragments_override", NULL)

if (is.null(rds) || is.null(out)) stop("Usage: --rds <obj.rds> --out <dir>")

if (!requireNamespace("Seurat", quietly=TRUE)) stop("Missing R package: Seurat")
if (!requireNamespace("Matrix", quietly=TRUE)) stop("Missing R package: Matrix")
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("Missing R package: jsonlite")
if (export_atac == "yes" && !requireNamespace("Signac", quietly=TRUE)) {
  stop("Missing R package: Signac")
}

dir.create(out, recursive=TRUE, showWarnings=FALSE)
obj <- readRDS(rds)

# group = active.ident
idents <- Seurat::Idents(obj)
barcodes <- names(idents)
groups <- as.character(idents)

write.table(data.frame(barcode=barcodes, group=groups),
            file=file.path(out,"cells.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)

# meta
meta <- obj@meta.data
meta$active_ident <- groups
con <- gzfile(file.path(out,"meta.tsv.gz"), "wb")
write.table(meta, file=con, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
close(con)

# reductions: prefer wnn.umap then umap
reductions <- Seurat::Reductions(obj)
pick_umap <- NULL
if ("wnn.umap" %in% reductions) pick_umap <- "wnn.umap"
if (is.null(pick_umap) && "umap" %in% reductions) pick_umap <- "umap"
if (!is.null(pick_umap)) {
  emb <- Seurat::Embeddings(obj, reduction=pick_umap)
  write.table(data.frame(barcode=rownames(emb), UMAP_1=emb[,1], UMAP_2=emb[,2]),
              file=file.path(out,"umap.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

assays <- Seurat::Assays(obj)
rna_assay <- if ("RNA" %in% assays) "RNA" else if ("SCT" %in% assays) "SCT" else NULL

# ---- TF expression summary (for gating) ----
# avg_expr from normalized data
if (!is.null(rna_assay)) {
  # avg_expr (data slot)
  avg_list <- Seurat::AverageExpression(obj, assays=rna_assay, slot="data", group.by="ident", verbose=FALSE)
  avg_mat <- avg_list[[rna_assay]]
  # pct_expr from raw counts
  cnt <- Seurat::GetAssayData(obj, assay=rna_assay, slot="counts")
  # group index
  grp_levels <- sort(unique(groups))
  pct_mat <- matrix(0, nrow=nrow(cnt), ncol=length(grp_levels))
  rownames(pct_mat) <- rownames(cnt)
  colnames(pct_mat) <- grp_levels
  for (g in grp_levels) {
    idx <- which(groups == g)
    if (length(idx) == 0) next
    sub <- cnt[, idx, drop=FALSE]
    pct_mat[, g] <- Matrix::rowSums(sub > 0) / length(idx)
  }

  # write long format: gene, group, avg_expr, pct_expr
  # (파일 사이즈 줄이기 위해 wide도 가능하지만, long이 파싱이 쉬움)
  genes <- rownames(avg_mat)
  long <- do.call(rbind, lapply(colnames(avg_mat), function(g) {
    data.frame(
      gene=genes,
      group=g,
      avg_expr=as.numeric(avg_mat[, g]),
      pct_expr=as.numeric(pct_mat[, g]),
      stringsAsFactors=FALSE
    )
  }))
  write.table(long, file=file.path(out,"tf_expr_by_group.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}

# ---- ATAC peaks + fragments map ----
if (export_atac == "yes") {
  atac_assay <- NULL
  if ("ATAC" %in% assays) atac_assay <- "ATAC"
  if (is.null(atac_assay) && "peaks" %in% assays) atac_assay <- "peaks"
  if (is.null(atac_assay)) stop("ATAC assay not found (expected ATAC or peaks).")

  feats <- rownames(obj[[atac_assay]])
  parse_one <- function(s) {
    if (grepl(":", s)) {
      parts <- strsplit(s, ":", fixed=TRUE)[[1]]
      chr <- parts[1]
      se <- strsplit(parts[2], "-", fixed=TRUE)[[1]]
      return(c(chr, se[1], se[2]))
    } else {
      se <- strsplit(s, "-", fixed=TRUE)[[1]]
      if (length(se) >= 3) return(c(se[1], se[2], se[3]))
    }
    return(c(NA, NA, NA))
  }
  mat <- do.call(rbind, lapply(feats, parse_one))
  colnames(mat) <- c("chr","start","end")
  mat <- mat[complete.cases(mat),,drop=FALSE]
  start0 <- as.integer(mat[, "start"]) - 1L
  start0[start0 < 0] <- 0
  bed <- data.frame(chr=mat[, "chr"], start=start0, end=as.integer(mat[, "end"]))
  write.table(bed, file=file.path(out,"peaks.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  frag_objs <- NULL
  try({ frag_objs <- Signac::Fragments(obj[[atac_assay]]) }, silent=TRUE)

  if (!is.null(frag_objs) && length(frag_objs) > 0) {
    mapping <- rep(NA_character_, length(barcodes)); names(mapping) <- barcodes
    for (f in frag_objs) {
      fp <- NA_character_; cells <- NULL
      try({ fp <- f@path }, silent=TRUE)
      try({ cells <- f@cells }, silent=TRUE)
      if (is.na(fp) || is.null(cells)) next
      cells <- intersect(cells, barcodes)
      mapping[cells] <- fp
    }
    if (any(is.na(mapping))) {
      if (!is.null(fragments_override)) mapping[is.na(mapping)] <- fragments_override
      else stop("Some cells missing fragment mapping. Provide --fragments_override.")
    }
    write.table(data.frame(barcode=names(mapping), fragment_path=as.character(mapping)),
                file=file.path(out,"fragments_map.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)
    writeLines(unique(as.character(mapping)), con=file.path(out,"fragments_unique.txt"))
  } else {
    if (!is.null(fragments_override)) {
      write.table(data.frame(barcode=barcodes, fragment_path=rep(fragments_override, length(barcodes))),
                  file=file.path(out,"fragments_map.tsv"),
                  sep="\t", quote=FALSE, row.names=FALSE)
      writeLines(fragments_override, con=file.path(out,"fragments_unique.txt"))
    } else {
      stop("No fragments found in object. Provide --fragments_override.")
    }
  }
}

# RNA counts export (for SCENIC+ bootstrap)
if (export_rna == "yes") {
  if (is.null(rna_assay)) stop("RNA assay not found (expected RNA or SCT).")
  mat <- Seurat::GetAssayData(obj, assay=rna_assay, slot="counts")
  Matrix::writeMM(mat, file.path(out, "rna_counts.mtx"))
  writeLines(rownames(mat), con=file.path(out, "genes.tsv"))
  writeLines(colnames(mat), con=file.path(out, "barcodes.tsv"))
}

report <- list(rds=rds, genome=genome, assays=assays, reductions=reductions, picked_umap=pick_umap,
               outputs=list(cells="cells.tsv", meta="meta.tsv.gz", umap="umap.tsv",
                            peaks="peaks.bed", fragments="fragments_map.tsv",
                            tf_expr="tf_expr_by_group.tsv"))
jsonlite::write_json(report, file.path(out,"export_report.json"), pretty=TRUE, auto_unbox=TRUE)
cat("[OK] export complete:", out, "\n")
