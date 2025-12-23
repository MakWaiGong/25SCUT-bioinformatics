# --------------- 环境初始化与依赖包管理 --------------- #
# 清空当前R环境所有变量，避免残留数据干扰
rm(list=ls())

# 验证Rtools是否安装（编译Bioconductor包必需）
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::find_rtools()  # 输出TRUE则Rtools安装成功

# 定义所需包列表（CRAN基础包 + Bioconductor生信分析包）
required_packages <- c(
  # CRAN包：数据处理/文件解压
  "tidyverse", "R.utils",
  # Bioconductor包：基因组分析核心工具
  "GenomicRanges", "IRanges", "rtracklayer", "biomaRt",
  # ChIP-seq分析专用包：峰注释/富集分析/差异分析
  "ChIPseeker", "clusterProfiler", "limma", 
  # 人类hg19参考基因组注释包（基因位置/ID映射）
  "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"
)

# 批量安装/加载包函数（自动判断CRAN/Bioconductor来源）
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # 区分Bioconductor包和CRAN包
      bioc_pkgs <- c("GenomicRanges", "IRanges", "rtracklayer", "biomaRt",
                     "ChIPseeker", "clusterProfiler", "limma", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
      if (pkg %in% bioc_pkgs) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(pkg, update = TRUE, ask = FALSE)  # 安装Bioconductor包
      } else {
        install.packages(pkg, dependencies = TRUE, ask = FALSE)  # 安装CRAN包
      }
    }
    library(pkg, character.only = TRUE)  # 加载包
  }
}
install_if_missing(required_packages)  # 执行包安装/加载

# --------------- 路径配置 --------------- #
# 原始数据文件夹（存储压缩的peaks.bed文件）
data_path <- "D:/document/code/lab7/Chip-seq/data"
# 结果输出文件夹（自动创建，存储分析结果）
output_path <- "D:/document/code/lab7/Chip-seq/results"
if (!dir.exists(output_path)) dir.create(output_path)

# --------------- 批量解压压缩文件 --------------- #
# 筛选所有peaks.bed.txt.gz格式的压缩文件
gz_files <- list.files(data_path, pattern = "peaks.bed.txt.gz", full.names = TRUE)
# 循环解压（仅当解压文件不存在时执行，避免重复操作）
for (gz_file in gz_files) {
  unzip_file <- gsub(".gz", "", gz_file)  # 解压后文件名（去除.gz后缀）
  if (!file.exists(unzip_file)) {
    gunzip(gz_file, destname = unzip_file, overwrite = FALSE, remove = FALSE)
    cat(paste("已解压：", basename(gz_file), "\n"))
  } else {
    cat(paste("已存在解压文件：", basename(unzip_file), "\n"))
  }
}

# --------------- 导入配对样本AR峰数据 --------------- #
# 定义3对配对样本的文件映射（正常组织vs肿瘤组织，一一对应）
sample_pairs <- list(
  patient1 = list(
    normal = "GSM1358395_DF_1335_normal_peaks.bed.txt",
    tumor = "GSM1358396_DF_1335_tumor_peaks.bed.txt"
  ),
  patient2 = list(
    normal = "GSM1358397_DF_1345_normal_peaks.bed.txt",
    tumor = "GSM1358398_DF_1345_tumor_peaks.bed.txt"
  ),
  patient3 = list(
    normal = "GSM1358399_DF_1373_normal_peaks.bed.txt",
    tumor = "GSM1358400_DF_1373_tumor_peaks.bed.txt"
  ),
  patient4 = list(
    normal = "GSM1358402_DF_1433_normal_peaks.bed.txt",
    tumor = "GSM1358403_DF_1433_tumor_peaks.bed.txt"
  ),
  patient5 = list(
    normal = "GSM1358405_DF_1609_normal_peaks.bed.txt",
    tumor = "GSM1358406_DF_1609_tumor_peaks.bed.txt"
  ),
  patient6 = list(
    normal = "GSM1358409_DF_184_normal_peaks.bed.txt",
    tumor = "GSM1358410_DF_184_tumor_peaks.bed.txt"
  )
)

# 批量导入bed格式的AR峰数据（存储为GenomicRanges对象，便于基因组区间分析）
ar_peaks <- list()
for (pair_name in names(sample_pairs)) {
  # 导入正常组织AR峰
  normal_file <- file.path(data_path, sample_pairs[[pair_name]]$normal)
  ar_peaks[[pair_name]]$normal <- import.bed(normal_file)
  # 导入肿瘤组织AR峰
  tumor_file <- file.path(data_path, sample_pairs[[pair_name]]$tumor)
  ar_peaks[[pair_name]]$tumor <- import.bed(tumor_file)
  # 打印导入信息（验证峰数量，确保导入成功）
  cat(paste(pair_name, "导入完成：正常组织峰数=", length(ar_peaks[[pair_name]]$normal), 
            "，肿瘤组织峰数=", length(ar_peaks[[pair_name]]$tumor), "\n"))
}

# --------------- 获取目标基因hg19基因组区间 --------------- #
# 目标基因列表（需分析AR结合的基因）
target_genes <- c("KLK4", "mTOR", "P4HB", "PDIA5")

# 1. 连接hg19版本的Ensembl数据库（GRCh37，匹配数据的基因组版本）
mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "grch37.ensembl.org"  # 关键：指定hg19对应的数据库
)

# 2. 获取基因的TSS（转录起始位点）和染色体信息
gene_info <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "transcription_start_site"),
  filters = "hgnc_symbol",  # 按基因名筛选
  values = target_genes,
  mart = mart
)

# 3. 基因去重（保留每个基因的第一个TSS，不影响区间分析）
gene_info_unique <- gene_info %>%
  group_by(hgnc_symbol) %>%
  slice(1) %>%
  ungroup()

# 4. 构建TSS上下游50kb区间（ChIP-seq分析中基因调控区的常用范围）
target_gene_regions <- GRanges(
  seqnames = paste0("chr", gene_info_unique$chromosome_name),  # 补充chr前缀，匹配bed文件格式
  ranges = IRanges(
    start = gene_info_unique$transcription_start_site - 50000,  # 上游50kb
    end = gene_info_unique$transcription_start_site + 50000,    # 下游50kb
    names = gene_info_unique$hgnc_symbol  # 用基因名命名区间
  )
)

# 验证目标基因区间（确保坐标正确）
cat("目标基因TSS上下游50kb区间（hg19）：\n")
print(target_gene_regions)

# --------------- 核心分析1：AR峰-基因区间重叠分析 --------------- #
# 功能：统计每对样本中，正常/肿瘤组织的AR峰与目标基因区间的重叠数量（判断是否存在AR结合）
peak_gene_overlap <- function(ar_peaks, target_regions) {
  overlap_results <- list()
  for (pair_id in names(ar_peaks)) {
    # 正常组织：计算AR峰与基因区间的重叠
    normal_overlap <- findOverlaps(ar_peaks[[pair_id]]$normal, target_regions)
    normal_result <- data.frame(
      pair_id = pair_id,
      tissue_type = "normal",
      gene = names(target_regions)[subjectHits(normal_overlap)],  # 重叠的基因名
      ar_peak_count = length(normal_overlap)  # 重叠峰数量
    ) %>%
      group_by(pair_id, tissue_type, gene) %>%
      summarise(ar_peak_count = n(), .groups = "drop")  # 按基因统计峰数
    
    # 肿瘤组织：同正常组织逻辑
    tumor_overlap <- findOverlaps(ar_peaks[[pair_id]]$tumor, target_regions)
    tumor_result <- data.frame(
      pair_id = pair_id,
      tissue_type = "tumor",
      gene = names(target_regions)[subjectHits(tumor_overlap)],
      ar_peak_count = length(tumor_overlap)
    ) %>%
      group_by(pair_id, tissue_type, gene) %>%
      summarise(ar_peak_count = n(), .groups = "drop")
    
    # 合并同一患者的正常/肿瘤结果
    pair_result <- rbind(normal_result, tumor_result)
    overlap_results[[pair_id]] <- pair_result
  }
  do.call(rbind, overlap_results)  # 合并所有患者结果
}

# 执行重叠分析
overlap_summary <- peak_gene_overlap(ar_peaks, target_gene_regions)

# 补充无AR峰的基因（填充0，保证数据完整性）
full_overlap <- expand.grid(
  pair_id = names(ar_peaks),
  tissue_type = c("normal", "tumor"),
  gene = target_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(overlap_summary, by = c("pair_id", "tissue_type", "gene")) %>%
  mutate(ar_peak_count = replace(ar_peak_count, is.na(ar_peak_count), 0))

# 打印重叠分析结果
cat("AR峰-基因重叠分析汇总（峰数量）：\n")
print(full_overlap)

# --------------- 核心分析2：提取AR峰强度（量化结合程度） --------------- #
# 功能：提取目标基因区间内AR峰的强度（BED第5列，RPM标准化值），计算平均/总强度
extract_peak_intensity <- function(ar_peaks, target_regions) {
  intensity_results <- list()
  for (pair_id in names(ar_peaks)) {
    # 正常组织：提取重叠峰的强度
    normal_peaks <- ar_peaks[[pair_id]]$normal
    normal_overlap_idx <- queryHits(findOverlaps(normal_peaks, target_regions))
    normal_intensity <- if (length(normal_overlap_idx) > 0) {
      normal_peaks[normal_overlap_idx] %>%
        as.data.frame() %>%
        mutate(
          pair_id = pair_id,
          tissue_type = "normal",
          gene = names(target_regions)[subjectHits(findOverlaps(normal_peaks[normal_overlap_idx], target_regions))],
          peak_intensity = score  # 提取峰强度（RPM）
        ) %>%
        group_by(pair_id, tissue_type, gene) %>%
        summarise(
          avg_peak_intensity = mean(peak_intensity),  # 平均强度（反映结合强度）
          total_peak_intensity = sum(peak_intensity),  # 总强度（反映结合总量）
          .groups = "drop"
        )
    } else {
      # 无重叠峰时填充0
      data.frame(
        pair_id = pair_id,
        tissue_type = "normal",
        gene = target_genes,
        avg_peak_intensity = 0,
        total_peak_intensity = 0,
        stringsAsFactors = FALSE
      )
    }
    
    # 肿瘤组织：同正常组织逻辑
    tumor_peaks <- ar_peaks[[pair_id]]$tumor
    tumor_overlap_idx <- queryHits(findOverlaps(tumor_peaks, target_regions))
    tumor_intensity <- if (length(tumor_overlap_idx) > 0) {
      tumor_peaks[tumor_overlap_idx] %>%
        as.data.frame() %>%
        mutate(
          pair_id = pair_id,
          tissue_type = "tumor",
          gene = names(target_regions)[subjectHits(findOverlaps(tumor_peaks[tumor_overlap_idx], target_regions))],
          peak_intensity = score
        ) %>%
        group_by(pair_id, tissue_type, gene) %>%
        summarise(
          avg_peak_intensity = mean(peak_intensity),
          total_peak_intensity = sum(peak_intensity),
          .groups = "drop"
        )
    } else {
      data.frame(
        pair_id = pair_id,
        tissue_type = "tumor",
        gene = target_genes,
        avg_peak_intensity = 0,
        total_peak_intensity = 0,
        stringsAsFactors = FALSE
      )
    }
    
    # 合并同一患者的强度结果
    pair_intensity <- rbind(normal_intensity, tumor_intensity)
    intensity_results[[pair_id]] <- pair_intensity
  }
  do.call(rbind, intensity_results)  # 合并所有患者结果
}

# 执行峰强度提取
intensity_summary <- extract_peak_intensity(ar_peaks, target_gene_regions)

# 合并峰数量+强度数据（完整AR结合信息）
combined_analysis <- full_overlap %>%
  left_join(intensity_summary[, c("pair_id", "tissue_type", "gene", "avg_peak_intensity", "total_peak_intensity")],
            by = c("pair_id", "tissue_type", "gene")) %>%
  mutate(
    avg_peak_intensity = replace(avg_peak_intensity, is.na(avg_peak_intensity), 0),
    total_peak_intensity = replace(total_peak_intensity, is.na(total_peak_intensity), 0)
  )

# 打印完整结合信息
cat("\nAR结合完整分析结果（峰数量+强度）：\n")
print(combined_analysis)

# --------------- 核心分析3：配对样本差异分析（肿瘤vs正常） --------------- #
# 预处理：计算肿瘤vs正常的强度倍数变化（FC），+0.1避免除0
diff_analysis_long <- combined_analysis %>%
  mutate(
    sample_id = paste(pair_id, tissue_type, sep = "_"),
    fc_avg_intensity = ifelse(tissue_type == "tumor",
                              (avg_peak_intensity + 0.1) / (avg_peak_intensity[match(paste(pair_id, "normal", sep = "_"), sample_id)] + 0.1),
                              NA)
  ) %>%
  filter(!is.na(fc_avg_intensity) | tissue_type == "normal")

# 按基因进行配对差异检验（limma包，FDR校正）
diff_results <- list()
for (gene in target_genes) {
  # 提取单个基因的所有样本数据
  gene_data <- diff_analysis_long %>%
    filter(gene == !!gene) %>%
    arrange(pair_id, tissue_type)
  
  # 无变异数据（全0）：直接生成无差异结果
  if (all(gene_data$avg_peak_intensity == 0)) {
    gene_diff <- data.frame(
      gene = gene,
      avg_fc = mean(gene_data$fc_avg_intensity[gene_data$tissue_type == "tumor"], na.rm = TRUE),
      logFC = 0,  # 无差异
      adj.P.Val = 1,  # 无显著差异
      total_pair_count = length(unique(gene_data$pair_id)),
      tumor_peak_pair = sum(gene_data$tissue_type == "tumor" & gene_data$ar_peak_count > 0),  # 肿瘤有峰的患者数
      normal_peak_pair = sum(gene_data$tissue_type == "normal" & gene_data$ar_peak_count > 0),  # 正常有峰的患者数
      peak_detail = paste(paste(gene_data$pair_id, gene_data$tissue_type, gene_data$ar_peak_count, sep = ":"), collapse = ",")
    )
    diff_results[[gene]] <- gene_diff
    next
  }
  
  # 构建强度矩阵（适配limma输入格式）
  intensity_matrix <- t(as.matrix(gene_data$avg_peak_intensity))
  rownames(intensity_matrix) <- "avg_peak_intensity"
  colnames(intensity_matrix) <- gene_data$sample_id
  
  # 定义配对信息和组织类型
  pair_info <- factor(gene_data$pair_id)
  tissue_info <- factor(gene_data$tissue_type, levels = c("normal", "tumor"))
  
  # 构建设计矩阵（肿瘤vs正常的对比）
  design <- model.matrix(~ tissue_info)
  colnames(design) <- c("Intercept", "Tumor_vs_Normal")
  
  # 校正配对样本的相关性
  corfit <- duplicateCorrelation(intensity_matrix, design, block = pair_info)
  if (is.na(corfit$consensus)) corfit$consensus <- 0  # 异常值保护
  
  # 线性模型拟合+Empirical Bayes校正（ChIP-seq差异分析标准方法）
  fit <- lmFit(intensity_matrix, design, block = pair_info, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  # 提取差异结果（FDR校正P值）
  gene_diff <- topTable(fit, coef = "Tumor_vs_Normal", number = Inf, adjust = "fdr") %>%
    mutate(
      gene = gene,
      avg_fc = mean(gene_data$fc_avg_intensity[gene_data$tissue_type == "tumor"], na.rm = TRUE),
      total_pair_count = length(unique(gene_data$pair_id)),
      tumor_peak_pair = sum(gene_data$tissue_type == "tumor" & gene_data$ar_peak_count > 0),
      normal_peak_pair = sum(gene_data$tissue_type == "normal" & gene_data$ar_peak_count > 0),
      peak_detail = paste(paste(gene_data$pair_id, gene_data$tissue_type, gene_data$ar_peak_count, sep = ":"), collapse = ",")
    ) %>%
    dplyr::select(gene, avg_fc, logFC, adj.P.Val, total_pair_count, tumor_peak_pair, normal_peak_pair, peak_detail)
  
  diff_results[[gene]] <- gene_diff
}

# 合并所有基因的差异结果
final_diff <- do.call(rbind, diff_results)

# 打印差异分析结果
cat("\n配对样本差异检验结果（肿瘤vs正常）：\n")
print(final_diff)

# --------------- 核心结论：AR调控判定规则 --------------- #
# 1. 有AR结合证据：基因TSS上下游50kb内存在AR峰（tumor_peak_pair≥1 或 normal_peak_pair≥1）；
# 2. 差异显著：FDR校正P值（adj.P.Val）≤0.05（小样本适配阈值，文献通常≤0.001）；
# 3. 肿瘤特异性结合（T-ARBS）：avg_fc>1 + tumor_peak_pair≥2 → 可能正调控；
# 4. 正常特异性结合（N-ARBS）：avg_fc<1 + normal_peak_pair≥2 → 可能负调控；
# 5. 无峰/无差异：无明确AR调控证据。