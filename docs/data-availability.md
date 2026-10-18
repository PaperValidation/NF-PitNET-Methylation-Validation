# Data Availability & Data-Source Notes

> 数据来源与说明的权威文档。全部原始数据为第三方公开数据，本仓库**不重分发原始数据**（仅代码 + 处理后的结果）。
> 对应稿件：`docs/manuscript-v2.md` 的 Data Availability 章节。

## 0. 总述

本研究是纯 **in-silico 二次分析**：所有原始甲基化/转录组数据来自公共数据库，研究者未生成任何新的人体样本数据。六个数据来源覆盖：

- **4 个甲基化队列**（3× EPIC-850K + 1× HM450）：用于跨队列、跨平台验证 NF-PitNET → GSE283928 c0（SF1-lineage）的映射。
- **1 个甲基化队列的原始 idat + 随访标签**（Møller 2025，Mendeley）：用于零 k-label 的从头投射与直接预后验证（再手术/再生长硬终点）。
- **1 个转录组队列**（GSE209903）：用于甲基化–转录一致性交叉验证。

### 本地存储约定（2026-10-18 更新，见项目 `.env`）

> **⚠️ 磁盘已于 2026-09~10 重组**：原 `G:\CMCC Sync\...` 目录已完全不存在，请勿再使用旧路径。

| 变量 | 路径 | 内容 | 备注 |
|---|---|---|---|
| `DATASET_PATH` | `E:\Medical Data\Neurosurgery\NF-PitNET` | **原始数据**：series matrix、idat、补充包 | `E:\Medical Data` 按科室分目录 |
| `DATASET_PATH_TEMP` | `F:\Data Temp\NF-PitNET` | **派生/中间/临时**：results、metadata、samples、intermediates、logs、papers、_wsl_scripts | 本地盘，不同步 |
| 项目 `results/` | `C:\...\Neurosurgery-NF-PitNET\results` | 投稿交付物（图、验证表、汇总 JSON） | — |

WSL 形态：`/mnt/e/Medical Data/Neurosurgery/NF-PitNET`、`/mnt/f/Data Temp/NF-PitNET`。

所有脚本通过 `scripts/_paths.py`（Python）/ `scripts/paths.R`（R）从 `.env` 读取路径，不硬编码。

---

## 1. 队列明细

### 1.1 Neou 2020 — E-MTAB-7762（EPIC-850K，训练/内部复现）

- **出处**：ArrayExpress/BioStudies — <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7762>
- **文献**：Neou M, et al. Pituitary neuroendocrine tumors. *Nat Commun* 2020.
- **平台**：Illumina Infinium MethylationEPIC 850K（**866,238 探针**；SDRF 标注的 450K 是错误标注——本分析按 EPIC 处理）。
- **样本**：86 例（32 NF-PitNET + 54 functioning）。
- **原始文件**：`.idat`（`idat_flat/`）；本地 `DATASET_PATH/E-MTAB-7762/`。
- **本研究的角色**：独立 EPIC 复现（81.2% → c0）；四种归一化敏感性分析（Table 3）；内部三分结构聚类（NFPT/TPIT/PIT1）。
- **注意**：EPIC idat 用 `minfi::read.metharray(force=TRUE)` 处理；本地派生 β 矩阵 `beta_out/neou_beta_cg.rds`（可再生成，不入库）。

### 1.2 GSE283928（Belakhoua 2025）— EPIC-850K（参考质心来源）

- **出处**：GEO — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283928>
- **文献**：Belakhoua S, et al. DNA methylation profiling of pituitary neuroendocrine tumors. *Neuro-Oncology* 2025;27(9):2341-2354. PMID 40295206.
- **平台**：EPIC-850K（submitter 处理后的 β 矩阵）。
- **样本**：118 例 PitNET 肿瘤；series matrix 含 **118 个真实 β 列 + 118 个 detection-pvalue 列（伪影，须剔除）**。剔除后无监督聚类得到 4 个肿瘤簇：**c0 n=70, c2 n=2, c3 n=22, c4 n=24**。
- **原始文件**：`GSE283928_MatrixProcessed_Normalized.txt.gz`（≈917 MB，本地 `DATASET_PATH/` 根目录）。
- **本研究的角色**：提供最近质心目标（c0/c2/c3/c4），全部下游投射都以它为参照。
- **注意**：原 milestone 中"cluster 1 = 正常垂体对照"是 detection-pvalue 伪影导致的错误结论，已订正（见 manuscript §3.1 与 Limitations）。

### 1.3 GSE207937（Silva-Junior 2023）— EPIC-850K（外部验证）

- **出处**：GEO — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207937>
- **文献**：Silva-Junior R, et al.（2023）。
- **平台**：EPIC-850K（submitter 处理后的 β 矩阵）。
- **样本**：77 例（45 非分泌型 NS/NF-PitNET + 32 functioning）。
- **原始文件**：`GSE207937_series_matrix.txt.gz`（本地 `DATASET_PATH/` 根目录）。
- **本研究的角色**：独立 EPIC 验证（75.6% → c0）。

### 1.4 GSE115783（Kober/Bujko 2018/2019）— HM450（跨平台外部验证）

- **出处**：GEO — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115783>
- **平台**：Illumina HumanMethylation450（HM450）——本研究的**跨平台**验证队列。
- **样本**：40 例 = **34 NFPA 肿瘤 + 6 正常垂体对照**（内嵌特异性对照，证明 c0 是肿瘤身份而非泛垂体信号）。
- **原始文件**：`GSE115783_series_matrix.txt.gz`（series matrix 即 ChAMP 处理后的 **β 值**，无需解 IDAT）；本地 `DATASET_PATH/GSE115783/`。
- **本研究的角色**：跨平台验证（94.1% → c0；正常对照仅 33.3% → c0）。

### 1.5 Møller 2025 — Mendeley 原始 idat（从头验证 + 直接预后）

- **出处**：Mendeley Data — DOI **10.17632/vr2wb38r6t.1** — <https://data.mendeley.com/datasets/vr2wb38r6t/1>
- **文献**：Møller J, et al.（2025）——提出五簇甲基化分型（k1–k5）与最高风险 k4 的原始队列。
- **平台**：EPIC-850K。
- **样本**：117 例 silent gonadotroph NF-PitNET（k1–k5 标签 + **纵向随访**：`ReOP`/`Reintervention` 与 `Cumulated regrowth %`）。
- **原始文件**：`Deidentified_IDATs/`（192 阵列 / 384 个 `.idat`）+ `SampleSheet.xlsx`；本地 `DATASET_PATH/Mendeley_Moller2025/`。
- **本研究的角色**：
  1. **零 k-label 从头验证**：用我们自己的 GSE283928 c0/c2/c3/c4 质心做最近质心投射（101/117 = 86.3% → c0）。
  2. **直接预后验证**：用其自带 k 标签与再手术硬终点（k3/k4/k5 vs k1/k2：HR 2.14, p=0.072；校正后 HR 2.09, p=0.088）。
  3. **562-DMP 探针级验证**：562 个已发表 cluster-specific DMP 套到独立 minfi 处理（99.1% LOO 恢复 k 分级）。
- **补充材料**：`noaf269_supplementary_data.zip`（PMC13070495 补充包）→ 本地 `DATASET_PATH/Moller2025_supplement/`，含 562-DMP 表（S2）与队列随访表（S1）。

### 1.6 GSE209903（Silva-Junior 2023）— 转录组（交叉验证）

- **出处**：GEO — <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209903>
- **平台**：RNA 转录组（差异表达）。
- **文件**：5 个对比表（FPT vs NFPT、A1 vs B1、A1 vs C1、B1 vs C1）+ series matrix；本地 `DATASET_PATH/` 根目录。
- **本研究的角色**：Møller 关键基因在 NF-PitNET 中启动子高甲基化 ↔ 转录下调的一致性检查（transcript_methyl_cross）。

---

## 2. 文件类型说明

| 类型 | 说明 | 处理方式 |
|---|---|---|
| `.idat`（EPIC/HM450） | 芯片原始荧光强度 | `minfi::read.metharray(force=TRUE)` → β |
| `*_series_matrix.txt.gz` | GEO submitter 处理矩阵 | 直接读（注意：GSE115783 是 β 值；GSE283928 需剔除 detection-pvalue 列） |
| `Mendeley_*` | Møller 原始 idat + SampleSheet | 同上（minfi force=TRUE，117 例 EPIC） |
| `noaf269_supplementary_data.zip` | Møller 补充表（562-DMP、随访） | 解压后读 xlsx/csv |

## 3. 引用语句（可直接用于稿件 Data Availability）

> All raw and processed methylation data analyzed in this study are publicly available. EPIC-850K datasets were obtained from GEO: GSE283928 (Belakhoua et al., 2025), GSE207937 (Silva-Junior et al., 2023), and GSE115783 (Kober/Bujko et al., 2018/2019); and from ArrayExpress/BioStudies: E-MTAB-7762 (Neou et al., 2020). Moller et al.'s raw EPIC 850K idat files were obtained from Mendeley Data (DOI 10.17632/vr2wb38r6t.1). Transcriptomic data were obtained from GSE209903. All custom analysis scripts and processed results are available in the project repository (github.com/tiduslrg/NF-PitNET-Methylation-Validation) and archived under Zenodo DOI 10.5281/zenodo.21840024 (restricted until acceptance).
