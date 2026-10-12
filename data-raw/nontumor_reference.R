# Builds inst/extdata/nontumor_genes.csv for screenNonTumor() (2026-10-11).
#
# Source: Human Protein Atlas, single-cell type data (rna_single_cell_type.tsv.zip, nCPM per gene and cell type;
# https://www.proteinatlas.org/about/download, CC BY-SA 4.0), downloaded 2026-10-11.
#
# A gene is listed for a non-tumor group if, in at least one of the group's HPA cell types, it has nCPM >= 10 and is
# >= 5x its highest nCPM in any epithelial or other tumor-like HPA cell type (so genes that carcinoma cells commonly
# express are not listed). For a site group, the site's own parenchymal types are excluded from the tumor-like set.
# Cell-cycle genes (Hallmark E2F targets, G2M checkpoint, mitotic spindle; the package's bundled gmt) are never
# listed: proliferating progenitors would otherwise contribute them, and tumor cells proliferate too.
# Usage: Rscript data-raw/nontumor_reference.R <path to rna_single_cell_type.tsv.zip>
suppressMessages(library(data.table))
zip <- commandArgs(trailingOnly = TRUE)[1]
d <- fread(cmd = paste("unzip -p", shQuote(zip)))
setnames(d, c("ens", "gene", "cell", "ncpm"))
m <- dcast(d, gene ~ cell, value.var = "ncpm", fun.aggregate = max, fill = 0)
m <- as.matrix(m[, -1], rownames = m$gene); m[!is.finite(m)] <- 0

GROUPS <- list(
  immune      = c("t-cells", "nk-cells", "b-cells", "plasma cells", "macrophages", "monocytes", "neutrophils",
                  "cdc", "pdcs", "mast cells"),
  stroma      = c("fibroblasts", "smooth muscle cells", "vascular smooth muscle cells", "pericytes"),
  vascular    = c("vascular endothelial cells", "lymphatic endothelial cells"),
  liver       = c("hepatocytes"),
  lung        = c("alveolar cells type 1", "alveolar cells type 2"),
  pancreas    = c("pancreatic acinar cells", "pancreatic islet cells"),
  brain       = c("astrocytes", "oligodendrocytes", "oligodendrocyte progenitor cells", "brain excitatory neurons",
                  "brain inhibitory neurons", "microglia"),
  kidney      = c("proximal tubule cells", "podocytes", "loop of henle epithelial cells", "distal convoluted tubule cells"),
  adrenal     = c("adrenal cortex cells", "adrenal medulla cells"),
  bone_marrow = c("erythrocyte progenitors", "megakaryocyte progenitors", "neutrophil progenitors", "monocyte progenitors"))
TUMORLIKE <- grep(paste0("epith|duct|keratinocyte|basal|secretory|colonocyte|enterocyte|cholangiocyte|club|",
                         "ciliated|goblet|alveolar|acinar|glandular|urothelial|tubule|hepatocyte|foveolar|chief|",
                         "parietal|islet|trophoblast|tuft|paneth|ionocyte|luminal|hormone|lactating|myoepithelial|",
                         "mucous|neuroendocrine|amplifying|deuterosomal|principal|intercalated|podocyte|papillary|henle"),
                  colnames(m), value = TRUE)
stopifnot(all(unlist(GROUPS) %in% colnames(m)))
out <- do.call(rbind, lapply(names(GROUPS), function(g) {
  ct <- GROUPS[[g]]; ref <- apply(m[, setdiff(TUMORLIKE, ct), drop = FALSE], 1, max)
  top <- apply(m[, ct, drop = FALSE], 1, max)
  genes <- rownames(m)[top >= 10 & top >= 5 * pmax(ref, 1)]
  data.frame(gene = genes, group = g, max_ncpm = round(top[genes], 1), ratio_vs_tumorlike = round(top[genes] / pmax(ref[genes], 1), 1))
}))
gmt <- strsplit(readLines("inst/extdata/h.all.v2024.1.Hs.symbols.gmt.txt"), "\t")
cc <- unique(unlist(lapply(gmt[sapply(gmt, `[`, 1) %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                                                          "HALLMARK_MITOTIC_SPINDLE")], `[`, -(1:2))))
out <- out[!out$gene %in% cc, ]
out <- out[order(out$group, out$gene), ]
write.csv(out, "inst/extdata/nontumor_genes.csv", row.names = FALSE)
print(table(out$group))
