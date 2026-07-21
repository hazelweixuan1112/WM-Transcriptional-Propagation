setwd("/Volumes/KWX/Grad/Programs/2025/202505/20250529/activity_gene_analysis/res/10.GO_enrichment")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

# ----------------------------
# 0) 输入：你的两份文件（分别是 GCI+ / GCI-）
#    约定：第一列是 Symbol（如果你文件不是这样，告诉我列结构我立刻改）
# ----------------------------
gci_pos_path <- "/Volumes/KWX/Grad/Programs/2025/202505/20250529/activity_gene_analysis/res/05.new_gene_score_list/gene_score_activate_increase_coupling_kwx.csv"
gci_neg_path <- "/Volumes/KWX/Grad/Programs/2025/202505/20250529/activity_gene_analysis/res/05.new_gene_score_list/gene_score_activate_decrease_coupling_kwx.csv"

# ----------------------------
# 1) 读取基因列表（兼容：有无表头；只取第一列 Symbol）
# ----------------------------
read_gene_symbols <- function(fp){
  df <- tryCatch(read.csv(fp, header=TRUE, stringsAsFactors=FALSE),
                 error=function(e) read.csv(fp, header=FALSE, stringsAsFactors=FALSE))
  # 取第一列作为 SYMBOL
  sym <- df[[1]]
  sym <- as.character(sym)
  sym <- trimws(sym)
  sym <- sym[!is.na(sym) & sym != ""]
  unique(sym)
}

gci_pos_sym <- read_gene_symbols(gci_pos_path)
gci_neg_sym <- read_gene_symbols(gci_neg_path)

# ----------------------------
# 2) SYMBOL -> ENTREZID（KEGG/Reactome 更稳）
#    同时做一个 universe（控制背景基因集合）
#    论文里是 portal 输入基因集；你这里用“两个集合的并集”作为背景是最干净的控制方式之一
# ----------------------------
all_sym <- unique(c(gci_pos_sym, gci_neg_sym))

to_entrez <- function(symbols){
  m <- bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  unique(m$ENTREZID)
}

pos_entrez <- to_entrez(gci_pos_sym)
neg_entrez <- to_entrez(gci_neg_sym)
uni_entrez <- to_entrez(all_sym)

# ----------------------------
# 3) Pathway enrichment（只做 KEGG + Reactome）
#    多重校正：BH (= FDR-BH)
#    阈值：p.adjust < 0.05（用 pvalueCutoff 控制输出）
# ----------------------------
run_kegg <- function(gene, universe){
  enrichKEGG(
    gene         = gene,
    organism     = "hsa",
    keyType      = "kegg",
    universe     = universe,
    pAdjustMethod= "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
  )
}

run_reactome <- function(gene, universe){
  enrichPathway(
    gene         = gene,
    organism     = "human",
    universe     = universe,
    pAdjustMethod= "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable     = TRUE
  )
}

pos_kegg <- run_kegg(pos_entrez, uni_entrez)
neg_kegg <- run_kegg(neg_entrez, uni_entrez)

pos_reac <- run_reactome(pos_entrez, uni_entrez)
neg_reac <- run_reactome(neg_entrez, uni_entrez)

# ----------------------------
# 4) 导出 TSV
#    为了后面 Python 画 “GeneRatio / Count / p.adjust / Description”
#    我们把 KEGG 和 Reactome 合并，并加一个 Source 标记
# ----------------------------
as_df_with_source <- function(enrich_obj, source_name){
  df <- as.data.frame(enrich_obj)
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  
  df2 <- df[, c("Description", "GeneRatio", "Count", "p.adjust")]
  df2$Source <- source_name
  df2
}


pos_df <- rbind(
  as_df_with_source(pos_kegg, "KEGG"),
  as_df_with_source(pos_reac, "Reactome")
)

# neg_df <- rbind(
  # as_df_with_source(neg_kegg, "KEGG"),
  # as_df_with_source(neg_reac, "Reactome")
# )

write.table(pos_df, file="GCI_pos_pathway.tsv", sep="\t",
            quote=FALSE, row.names=FALSE)
# write.table(neg_df, file="GCI_neg_pathway.tsv", sep="\t",
            # quote=FALSE, row.names=FALSE)

# cat("Done!\nGenerated:\n  GCI_pos_pathway.tsv\n  GCI_neg_pathway.tsv\n")