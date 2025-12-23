# ============================================================================ #
# 蛋白质组分析####
# ============================================================================ #
# 环境准备####
#清空环境
rm(list=ls())
# 库,养成不超过细线的好习惯！
required_packages <- c(
  "tidyverse", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr",
  "DESeq2","tibble","limma","dplyr","tibble"
)

install_if_missing <- function(pkgs){
  for (pkg in pkgs){
    if (!requireNamespace(pkg, quietly = TRUE)){
      if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr")){
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg)
      } else install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

install_if_missing(required_packages)
# 显示声名随机种子，是保持复现性的好习惯
set.seed(42)

# 项目结构
# --
# |-proteomics_analysis.R
# |-transcriptomics_analysis.R
# |-data
# | |-Spectronaut_Output-Result.tsv 蛋白组数据
# | |-GSE288596_RPKM_allSamples.xlsx 转录组数据
# | |-hsa_hallmark_genes.rds
# |-result
# | |-transcriptomics
# | |-proteomics
setwd("D:/document/code/lab7")
output_path <- "result//proteomics"
dir.create("result", showWarnings = FALSE)
dir.create("result/proteomics", showWarnings = FALSE)
# ============================================================================ #
# 数据加载与解析####
# ============================================================================ #
raw_df <- read_tsv("data/Spectronaut_Output-Result.tsv", show_col_types = T)
# 查看raw_df以知道后续解析
#> raw_df
# # A tibble: 2,150 × 128
# PG.ProteinGroups     PG.Genes        PG.ProteinDescriptions  PG.ProteinNames PG.NrOfStrippedSeque…¹
# <chr>                <chr>           <chr>                   <chr>                            <dbl>
#   1 A0A0B4J2D5;P0DPI2 GATD3B;GATD3A   Glutamine amidotransfe… GAL3B_HUMAN;GA…                      5
# 2 A0AV96               RBM47           RNA-binding protein 47  RBM47_HUMAN                          1
# 3 A0FGR8               ESYT2           Extended synaptotagmin… ESYT2_HUMAN                          3
# 4 A0MZ66               SHTN1           Shootin-1               SHOT1_HUMAN                          1
# 5 A1L0T0               ILVBL           2-hydroxyacyl-CoA lyas… HACL2_HUMAN                          3

expr_df <- raw_df %>%
  mutate(Gene = str_split(PG.Genes, ";", simplify = TRUE)[, 1]) %>%
  # simplify = TRUE 返回字符矩阵而非列表便于操作
  dplyr::filter(!is.na(Gene)) %>% # 忘记dplyr::空值命名空间，是容易踩的坑
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  # .keep_all = TRUE保留其他行
  column_to_rownames("Gene")

expr_log2 <- expr_df[, grep("PG.Log2Quantity", colnames(expr_df), value = TRUE)]

# 分组函数
split_exp <- function(mat, pattern, label_fun){
  m <- mat[, grep(pattern, colnames(mat)), drop = FALSE]
  labels <- label_fun(colnames(m))
  list(mat = na.omit(m), labels = labels)
}
# 正则表达式CCF642|Vehicle的|是或的含义
sm_exp <- split_exp(expr_log2, "CCF642|Vehicle", 
                    function(x) ifelse(grepl("CCF642", x), "CCF642", "Vehicle"))

sirna_exp <- split_exp(expr_log2, "PDIA1|PDIA5|SC1",
                       function(x) case_when(
                         grepl("PDIA1", x) ~ "siPDIA1",
                         grepl("PDIA5", x) ~ "siPDIA5",
                         TRUE ~ "siNC"))
# ============================================================================ #
# 工具函数####
# ============================================================================ #
# SYMBOL -> ENTREZID
map_symbol_to_entrez <- function(symbols){
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}


# 差异分析 edgeR和DESeq2要count limma可直接使用log2 所以limma

diff_expression_limma <- function(mat, labels, g1, g2) {
  # 加载依赖包（确保函数内可调用）
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("请安装dplyr包：install.packages('dplyr')")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("请安装tibble包：install.packages('tibble')")
  if (!requireNamespace("limma", quietly = TRUE)) stop("请安装limma包：install.packages('limma')")
  
  # ========== 前置检验（核心：提前定位问题） ==========
  cat("=== 开始limma差异分析：", g1, "vs", g2, "===\n")
  
  # 检验1：初始mat列数和labels长度是否匹配
  if (ncol(mat) != length(labels)) {
    stop(paste0("初始数据不匹配！mat列数=", ncol(mat), "，labels长度=", length(labels)))
  }
  cat("  初始数据：mat列数=", ncol(mat), "，labels长度=", length(labels), "\n")
  
  # 检验2：g1/g2是否在labels中（避免分组名拼写错误）
  if (!g1 %in% labels | !g2 %in% labels) {
    stop(paste0("分组名错误！labels中无", ifelse(!g1%in%labels, g1, g2), "，请检查拼写"))
  }
  cat("  分组验证：", g1, "样本数=", sum(labels==g1), "，", g2, "样本数=", sum(labels==g2), "\n")
  
  # ========== 步骤1：正确筛选样本（按labels筛选，而非列名） ==========
  # 生成筛选索引：保留labels中属于g1/g2的样本
  keep_idx <- labels %in% c(g1, g2)
  # 同步筛选mat和labels（drop=FALSE避免单列降级）
  mat_filtered <- mat[, keep_idx, drop = FALSE]
  labels_filtered <- labels[keep_idx]
  
  # 检验3：筛选后维度
  cat("  筛选后：mat列数=", ncol(mat_filtered), "，labels长度=", length(labels_filtered), "\n")
  if (ncol(mat_filtered) == 0) stop("筛选后mat无列！请检查g1/g2是否正确")
  if (ncol(mat_filtered) != length(labels_filtered)) stop("筛选后维度不匹配")
  
  # ========== 步骤2：构建design矩阵（对照组g2在前，实验组g1在后） ==========
  # 因子化标签：参考水平=对照组g2
  group_factor <- factor(labels_filtered, levels = c(g2, g1))
  # 构建无截距设计矩阵
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- c(g2, g1)  # 列名对应对照组、实验组
  cat("  design矩阵：行数=", nrow(design), "，列数=", ncol(design), "\n")
  
  # ========== 步骤3：limma核心分析（补全contrasts.fit） ==========
  # 拟合线性模型
  fit <- limma::lmFit(mat_filtered, design)
  # 构建对比矩阵（实验组g1 - 对照组g2）
  contrast <- limma::makeContrasts(contrasts = paste0(g1, "-", g2), levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  # 贝叶斯方差收缩
  fit2 <- limma::eBayes(fit2)
  
  # ========== 步骤4：提取并整理结果 ==========
  # 提取所有基因的差异结果（coef=1对应g1-g2的对比列）
  res <- limma::topTable(fit2, coef = 1, number = Inf, adjust = "fdr")
  
  # 整理结果：行名转列、重命名、过滤无效值
  result_df <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%  # 基因名从行名转列
    dplyr::select(Gene, logFC, P.Value, adj.P.Val) %>%  # 保留原始P值和校正P值
    dplyr::rename(log2FC = logFC, pvalue = P.Value, padj = adj.P.Val) %>%  # 统一列名
    dplyr::filter(!is.na(Gene) & !is.na(padj))  # 过滤无效值
  
  # ========== 步骤5：打印分析结果 ==========
  sig1 <- sum(result_df$padj < 0.05)  # padj<0.05的差异基因
  sig2 <- sum(result_df$padj < 0.01 & abs(result_df$log2FC) > 0.5)  # 严格显著的基因
  cat("  差异分析完成！\n")
  cat("  差异基因总数（padj<0.05）：", sig1, "\n")
  cat("  显著差异基因（padj<0.01 & |log2FC|>0.5）：", sig2, "\n")
  
  # 返回tibble格式（适配后续分析）
  return(tibble::as_tibble(result_df))
}


# ============================================================================ #
# Hallmark下载与缓存####
# ============================================================================ #
# 缓存路径和文件名
cache_path <- "data"  # 可自定义
if (!dir.exists(cache_path)) dir.create(cache_path, recursive = TRUE)
cache_file <- file.path(cache_path, "hsa_hallmark_genes.rds")

# 下载或读取Hallmark基因集
if (!file.exists(cache_file)) {
  hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, gene_symbol)
  saveRDS(hallmark_genes, cache_file)
  message("Hallmark基因集已下载并缓存至：", cache_file)
} else {
  hallmark_genes <- readRDS(cache_file)
  message("已从缓存读取Hallmark基因集：", cache_file)
}

# 查看基本信息--->找到文献要富集的ROS通路名
# > head(hallmark_genes)           # 前6行
# # A tibble: 6 × 2
# gs_name               gene_symbol
# <chr>                 <chr>      
#   1 HALLMARK_ADIPOGENESIS ABCA1      
# 2 HALLMARK_ADIPOGENESIS ABCB8      
# 3 HALLMARK_ADIPOGENESIS ACAA2      
# 4 HALLMARK_ADIPOGENESIS ACADL      
# 5 HALLMARK_ADIPOGENESIS ACADM      
# 6 HALLMARK_ADIPOGENESIS ACADS      
# > unique(hallmark_genes$gs_name) 
# [1] "HALLMARK_ADIPOGENESIS"                     
# [2] "HALLMARK_ALLOGRAFT_REJECTION"              
# [3] "HALLMARK_ANDROGEN_RESPONSE"                
# [4] "HALLMARK_ANGIOGENESIS"  
# ......
# [42] "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  --->找到文献要富集的ROS通路名
# ......      
# [47] "HALLMARK_UV_RESPONSE_DN"                   
# [48] "HALLMARK_UV_RESPONSE_UP"                   
# [49] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"       
# [50] "HALLMARK_XENOBIOTIC_METABOLISM"            
# > nrow(hallmark_genes)           # 总基因数
# [1] 7333

# SYMBOL -> ENTREZID 并生成 Hallmark list
hallmark_entrez <- hallmark_genes %>%
  inner_join(map_symbol_to_entrez(unique(hallmark_genes$gene_symbol)),
             by = c("gene_symbol" = "SYMBOL")) %>%
  distinct(gs_name, ENTREZID, .keep_all = TRUE)

hallmark_list <- split(hallmark_entrez$ENTREZID, hallmark_entrez$gs_name)

# ============================================================================ #
# GSEA分析####
# ============================================================================ #
run_gsea_hallmark <- function(de_df, comparison_name, pathways,out_dir="result/proteomics"){
  ros_pathway_name <- "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
  de_df <- de_df %>% distinct(Gene, .keep_all = TRUE)
  
  mapped <- bitr(de_df$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  gsea_df <- de_df %>% inner_join(mapped, by = c("Gene"="SYMBOL")) %>% distinct(ENTREZID, .keep_all = TRUE)
  
  ranks <- sort(setNames(gsea_df$log2FC, gsea_df$ENTREZID), decreasing = TRUE)
  
  term2gene <- map_dfr(names(pathways), function(term) tibble(TERM=term, GENE=pathways[[term]]))
  
  fg <- clusterProfiler::GSEA(
    geneList = ranks,
    TERM2GENE = term2gene,
    TERM2NAME = term2gene[, c("TERM","TERM")],
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "fdr",
    verbose = FALSE
  )
  fg@result$Description[fg@result$ID == ros_pathway_name] <- "ROS Pathway"
  fg_df <- as.data.frame(fg) %>% arrange(p.adjust)
  write_csv(fg_df,file.path(out_dir,paste0(comparison_name, "_gsea_hallmark.csv")) )
  
  # Top15柱状图
  p <- ggplot(fg_df %>% head(15), aes(reorder(Description, NES), NES, fill=NES)) +
    geom_col() + coord_flip() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    theme_minimal() +
    labs(title=paste("GSEA Hallmark:", comparison_name), x="", y="NES") +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  print(p)
  ggsave(file.path(out_dir,paste0(comparison_name, "_gsea_hallmark.tiff")), p, width=8, height=5, dpi=300)
 
  # 检查ROS通路是否在GSEA结果中
  if (!ros_pathway_name %in% fg$ID) {
    warning(paste(
      "ROS通路 [", ros_pathway_name, "] 未在", comparison_name, "的GSEA结果中找到！",
      "\n请检查通路名是否匹配，或是否有足够的有效基因"
    ))
  } else {
    # 绘制ROS通路的GSEA富集图（含ES曲线、基因条码、统计表格）
    ros_plot <- gseaplot2(
      fg,  # clusterProfiler返回的gseaResult对象（原生适配）
      geneSetID = ros_pathway_name,  # 指定ROS通路
      pvalue_table = TRUE,  # 显示NES/padj/pvalue统计表格
      color = "#E74C3C",    # 富集曲线颜色（可自定义）
      # title = paste("GSEA Enrichment: ROS Pathway (", comparison_name, ")"),
      base_size = 12        # 字体大小
    )
  }
    print(ros_plot)
    # 保存ROS通路富集图
    ggsave(
      file.path(out_dir,paste0(comparison_name, "_ROS_pathway_gsea_enrichment.tiff")),
      ros_plot,
      width = 10,  # 宽度适配富集图
      height = 6,
      dpi = 300
    )
    return(fg)
      }

# ============================================================================ #
# ORA分析 (GO的BP)####
# ============================================================================ #
run_ora_go <- function(de_df, universe_symbols, comparison,
                       ont = "BP", log2fc = log2(1.5), p_cutoff = 0.05,
                       min_genes = 5, show_n = 8, out_dir = "result/proteomics") {
  # 1. 筛选差异基因（SYMBOL）
  deg_genes <- de_df %>%
    filter(!is.na(Gene), abs(log2FC) >= log2fc, pvalue < p_cutoff) %>%
    pull(Gene) %>% unique()
  
  # 2. SYMBOL转ENTREZID（精简辅助函数）
  sym2entrez <- \(x) bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db) %>% 
    pull(ENTREZID) %>% unique() %>% na.omit()
  deg_entrez <- sym2entrez(deg_genes)
  uni_entrez <- sym2entrez(universe_symbols)
  
  # 3. GO富集分析（仅保留核心判断）
  if (length(deg_entrez) >= min_genes) {
    ego <- enrichGO(gene = deg_entrez, 
                    # universe = uni_entrez,
                    OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = ont, pAdjustMethod = "fdr", readable = TRUE)
    
    # 4. 结果输出（仅当富集有结果时执行）
    if (!is.null(ego) && nrow(ego@result) > 0) {
      go_df <- as_tibble(ego@result) %>% arrange(desc(FoldEnrichment))
      dir.create(out_dir, showWarnings = FALSE)
      
      # 保存CSV
      write_csv(go_df, file.path(out_dir, paste0(comparison, "_ORA_GO_", ont, ".csv")))

      # 绘图+保存
      plot_df <- go_df %>% slice_head(n = show_n)
      # 绘图前添加，查看qvalue实际分布
      cat("当前top通路qvalue范围：", min(plot_df$qvalue), "~", max(plot_df$qvalue), "\n")
      cat("qvalue中位数：", median(plot_df$qvalue), "\n")
      p <- ggplot(plot_df, aes(FoldEnrichment, reorder(Description, FoldEnrichment))) +
        geom_point(aes(size = Count, color = qvalue)) +
        scale_color_gradient(low = "red", high = "blue") +
        scale_color_viridis_c(option = "plasma", trans = "log10")+
        theme_bw() +
        labs(title = paste("GO ORA", ont, comparison, sep = " | "),
             x = "Fold Enrichment", y = NULL, color = "q value")
      
      ggsave(file.path(out_dir, paste0(comparison, "_ORA_GO_", ont, ".tiff")),
             p, width = 8, height = max(5, 0.25 * nrow(plot_df)), dpi = 600)
      print(p)
      
      return(invisible(go_df))
    }
  }
}

# ============================================================================ #
# 执行分析####
# ============================================================================ #
res_pdia1_limma <- diff_expression_limma(sirna_exp$mat, sirna_exp$labels, "siPDIA5", "siNC")
universe_symbols <- rownames(expr_log2)
gsea_pdia1 <- run_gsea_hallmark(res_pdia1_limma, "siPDIA5_vs_siNC", hallmark_list)
ora_pdia1 <- run_ora_go(res_pdia1_limma, universe_symbols, "siPDIA5_vs_siNC")

res_pdia1_limma <- diff_expression_limma(sirna_exp$mat, sirna_exp$labels, "siPDIA1", "siNC")
universe_symbols <- rownames(expr_log2)
gsea_pdia1 <- run_gsea_hallmark(res_pdia1_limma, "siPDIA1_vs_siNC", hallmark_list)
ora_pdia1 <- run_ora_go(res_pdia1_limma, universe_symbols, "siPDIA1_vs_siNC")

res_pdia1_limma <- diff_expression_limma(sm_exp$mat, sm_exp$labels, "CCF642", "Vehicle")
universe_symbols <- rownames(expr_log2)
gsea_pdia1 <- run_gsea_hallmark(res_pdia1_limma, "CCF642_vs_Vehicle", hallmark_list)
ora_pdia1 <- run_ora_go(res_pdia1_limma, universe_symbols, "CCF642_vs_Vehicle")

