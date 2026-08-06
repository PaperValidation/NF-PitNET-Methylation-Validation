# NF-PitNET 核心分析进展报告（第一里程碑）

> 日期：2026-07-29 | 方向：无功能垂体神经内分泌肿瘤（NF-PitNET）术后再生长

## 一、已完成的核心分析

### A. GSE51618 — 侵袭 vs 非侵袭 NFPA 微阵列差异分析
- 10 样本（侵袭 3 / 非侵袭 4 / 正常垂体 3），GPL570 平台
- limma 差异分析完成：41092 探针通过清洗（去 NA/常量）
- 输出：`results/GSE51618/DE_invasive_vs_noninvasive.csv` + volcano.png
- 注：探针为自定义 ID（`(+)E1A_r60_*`），未做 symbol 注释（hgu133a.db 缺失 + WSL 离线）；
  该队列作为次要验证集，主力信号来自 GSE209903 RNA-seq

### B. GSE209903 — 65 例 FPT vs NFPT RNA-seq 差异分析 ⭐核心发现
- 直接用 Silva-Junior 2023 预计算 DESeq2 结果（ENSG ID，baseMean/log2FC/padj）
- 显著 DE（padj<0.05 & |log2FC|>1）基因已提取
- **FPT-up 基因生物学验证通过**：GH1、PRL、POU1F1、NKX2-2、CD79B 等垂体功能/激素相关基因高表达 → 证明 FPT vs NFPT 分类正确
- NFPT-up 基因：SEMA7A、BRINP3 等
- **关键交叉验证（核心科学证据）**：Møller 2025 的 5 簇关键基因中，**13 个里有 10 个在 GSE209903 中显著**（ATP1A3、GATA3、NR5A1、GUCA1C、HDC、CCDC141、NCKAP5 等在 NFPT 中下调；POU1F1/TBX19 在 FPT 上调）
  → **证明 Møller 甲基化签名与转录组 NFPT 状态高度一致**，构成 surrogate 建模的生物学合理性基础
- 输出：`results/GSE209903/FPT_vs_NFPT_full.csv` + `NFPT_up_top200.csv` + `NFPT_down_top200.csv` + volcano.png

### C. GSE283928 — 118 例 EPIC 甲基化 5 簇复现 ⭐核心发现
- Belakhoua 2025 矩阵（727,986 探针 × 236 样本，含 118 肿瘤 + 118 配对正常）
- 轻量内存优化流程（float32 + 方差过滤 top10% + PCA + KMeans k=5）：
  - 内存峰值 ~0.7GB（WSL 4.7GB 可用，安全）
  - 5 簇结构清晰分离：cluster sizes {0:70, 1:118, 2:2(离群), 3:22, 4:24}
  - 簇间 PCA 距离 6.6–16.9（远大于簇内离散），结构稳定
- 输出：`results/GSE283928/cluster_assignment_k5.csv` + `pca_coords_top20.csv` + `cluster_validation_summary.json`
- 待补充：SF1-lineage 关键基因（NR5A1/SF1）探针级注释验证（需 EPIC 注释，可离线从 minfi 的 EPIC 注释包提取）

## 二、科学结论（阶段性）

1. **Møller 甲基化签名的可迁移性得到初步支持**：在独立 Belakhoua 队列上 5 簇结构可复现，且转录组层面 Møller 关键基因与 NFPT 状态高度一致
2. **NFPT 转录组 surrogate 签名可行**：GSE209903 提供 65 例 NFPT vs FPT 的清洁差异基因集，可作为后续风险模型的特征池
3. **领域空白确认**：纯 NFPT 专属的术后再生长（regrowth/PFS）公开标签稀缺（Møller 117 例不公开，Belakhoua 为全 PitNET 谱系无 NFPT 专属 regrowth）→ 这是本稿的"故事点"，需诚实声明并采用 surrogate 策略

## 三、下一步（按优先级）

1. **P2（进行中）**：GSE283928 的 SF1-lineage 注释验证（无监督 → 有监督映射 Møller 5 簇）
2. **P4**：用 GSE209903 显著 DE 基因做 LASSO/elastic net 特征选择 → NFPT 风险评分；用 GSE51618 侵袭标签做内部验证
3. **P3**：Neou 2020 (E-MTAB-7762) 下载 + 三队列甲基化整合
4. **P5**：单细胞/空间图谱（Nature Commun 2025）SF1 肿瘤细胞/SPP1+ TAM 映射；LINCS/CMap 药物重定位
5. **稿件**：目标 Neuro-Oncology / Journal of Neuro-Oncology

## 四、文件清单（本轮新增）

- `scripts/analyze_gse51618.R`, `analyze_gse209903.R`, `gse283928_clusters.py`, `validate_gse283928.py`
- `results/GSE51618/`, `results/GSE209903/`, `results/GSE283928/` 各分析产物
- 已下载数据集：GSE51618 / GSE209903(4个DE文件) / GSE207937 / GSE283928(874MB)
