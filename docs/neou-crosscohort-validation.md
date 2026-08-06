# Neou 2020 (E-MTAB-7762) 跨队列甲基化验证 — 里程碑 7 报告

**日期**：2026-07-30
**队列**：Neou 2020, *E-MTAB-7762* (Neou *et al.*, *Acta Neuropathol* 2020) — 86 例 EPIC 850K 甲基化（经 minfi 从 raw idat 计算 β）
**目的**：将第三个独立 EPIC 队列整合进 Møller 2025 / Belakhoua 2025 甲基化分型框架，验证该框架对 NF-PitNET 生物学的捕获能力。

---

## 1. 关键事实修正

### 1.1 Neou 平台是 EPIC 850K，不是 SDRF 写的 450K
- SDRF 中标注平台为 `A-GEOD-21145`（Illumina HumanMethylation 450K），但实际 idat 经 minfi 读取后得到 **866,238 probes**，为 EPIC 850K 特征。
- 因此 Neou 是 **第三个独立 EPIC 队列**，不是此前笔记中提到的“450K↔EPIC 跨平台”。跨队列一致性由 3 个 EPIC 队列共同支撑。

### 1.2 GSE283928 本地矩阵结构修正
- 本地 `GSE283928_MatrixProcessed_Normalized.txt.gz` 含 **118 个真实肿瘤 β 列**（`Case1_Beta` … `Case118_Beta`）+ **118 个 detection-pvalue 列**（`Case1_Detection Pval` … `Case118_Detection Pval`，数值 0–0.01，78.9% 为 0）。
- 此前里程碑把 cluster 1（n=118）误判为“正常垂体对照”，实际是 **detection-pvalue 列聚成的伪影**。
- 118 个真实肿瘤 β 列被正确划分为 4 个肿瘤簇：**c0:70, c2:2, c3:22, c4:24**。该 4 簇结构用于后续投影。

---

## 2. 方法

1. **数据形态**
   - Neou: `neou_beta_cg.rds` → 866,238 probes × 86 samples（cg-keyed，EPIC hg19 注释）。
   - GSE283928: 仅保留 118 个 `_Beta` 列 → 727,986 probes × 118 samples。
   - GSE207937: 系列矩阵 M 值 → `β = 1/(1+2^(-M))` → 769,976 probes × 77 samples。

2. **缺失值处理**
   - 每个矩阵对 probe 做列均值填充；若整列缺失则填 0.5。

3. **公共探针与联合 PCA**
   - 三队列交集：**687,463 common probes**。
   - 保留方差 top 10% → **68,746 probes**；PCA(20)（sklearn）。
   - PC1 = 29.2%, PC2 = 9.7%, cum20 = 63.4%。

4. **跨队列质心投影**
   - 在 687,463 common-probe β 空间中计算 GSE283928 4 个肿瘤簇质心。
   - 将 Neou 与 GSE207937 样本按**欧氏距离**和**余弦相似度**分配到最近质心，与 GSE207937 原分析方法保持一致。

5. **谱系 / Møller-2025 关键基因评分**
   - 使用 `results/GSE283928/lineage_marker_probes.csv` 中的 SF1_lineage / TPIT_lineage / PIT1_lineage / Moller_key 探针集。
   - 谱系评分：`SF1 - (TPIT + PIT1) / 2`。
   - 置换检验（20,000 次）比较 NFPT vs 分泌型。

---

## 3. 核心结果

### 3.1 联合 PCA：Neou 与既有队列共享同一甲基化流形

| 队列 | n | PC1 均值 | PC2 均值 | PC1 std | PC2 std |
|---|---|---|---|---|---|
| Neou | 86 | +8.40 | +1.50 | 38.92 | 20.36 |
| GSE283928 | 118 | −4.39 | −1.21 | 23.68 | 15.57 |
| GSE207937 | 77 | −2.65 | +0.18 | 26.99 | 16.85 |

队列质心间距离：
- GSE283928 ↔ GSE207937：**2.23**（两个参考队列高度重叠）
- Neou ↔ GSE207937：**11.13**
- Neou ↔ GSE283928：**13.08**

→ Neou 的 PC1/PC2 标准差（38.9 / 20.4）显著大于队列质心偏移，说明 Neou 样本**广泛分布于既有队列的甲基化流形中**，并非形成孤立的批次岛。

### 3.2 跨队列质心投影：NFPT 稳定映射到 cluster 0

**Neou**
| 亚组 | c0 | c2 | c3 | c4 |
|---|---|---|---|---|
| NFPT (n=32) | **26 (81%)** | 0 | 2 (6%) | 4 (13%) |
| Functioning (n=54) | 2 (4%) | 0 | 15 (28%) | **37 (69%)** |

**GSE207937**（重跑一致性核对）
| 亚组 | c0 | c2 | c3 | c4 |
|---|---|---|---|---|
| NS (NFPT, n=45) | **35 (78%)** | 0 | 6 (13%) | 4 (9%) |
| GH/ACTH/NA (n=32) | 0 | 0 | 7 | 25 |

→ 两个独立 NFPT 队列中，**绝大多数 NFPT 样本投影到 GSE283928 tumor cluster 0**；而分泌型/功能性肿瘤则倾向 cluster 3/4。这是跨队列可迁移性的直接证据。

### 3.3 Møller-2025 关键基因甲基化：NFPT 显著高甲基化

| 队列/亚组 | mean β @ Moller-key genes |
|---|---|
| Neou NFPT | **0.687 ± 0.075** |
| Neou Functioning | 0.563 ± 0.096 |
| GSE207937 NS (NFPT) | **0.629 ± 0.058** |
| GSE207937 other | 0.597 ± 0.052 |

- Neou NFPT vs Functioning：p ≈ 0.00005（置换检验）
- GSE207937 NS vs other：p ≈ 0.04（此前报告）

→ **Møller-2025 关键基因在 NFPT 中呈高甲基化，在两个独立队列中复现**；结合 GSE209903 转录组结果（这些基因在 NFPT 中转录下调），形成 **甲基化↑ ↔ 表达↓** 的跨组学一致证据。

### 3.4 谱系评分：方向一致但需谨慎解读

| 队列/亚组 | Lineage score (SF1 − (TPIT+PIT1)/2) |
|---|---|
| Neou NFPT | **−0.117** |
| Neou Functioning | −0.010 |
| GSE207937 NS (NFPT) | **−0.046** |
| GSE207937 other | +0.014 |

→ NFPT 在两个队列中均略低于分泌型，方向与“NFPT = SF1-lineage 主导”的简单预期不一致。这与此前 GSE207937 分析的结论相同：单纯用平均 β 差值构建的谱系代理评分过于粗糙，**不应作为独立生物学证据**。更稳健的锚点是 Møller 关键基因甲基化和跨队列投影。

### 3.5 Neou 内部独立聚类：从零复现跨队列结构（方案 C）

为验证上述结构是否依赖 GSE283928 质心、还是可由 Neou 自身数据独立产生，在 Neou 86 例 EPIC 样本上单独执行无监督聚类。

- 特征选择：Neou 内部方差 top 10% probes → **86,623 probes**。
- 降维：PCA(20), random_state=42；PC1=45.6%, PC2=8.6%, cum20=77.0%。
- 聚类：KMeans k=2..6，以 **轮廓系数** 选最优 k；k=3 最优（silhouette 0.4045）。
- 簇标签按 Møller-key 甲基化均值降序重排，使 cluster 0 = 最高风险样态。

| 内部簇 | n | NFPT 占比 | Møller-key β | 谱系评分 | 主要病理 | 跨队列投射 |
|---|---|---|---|---|---|---|
| **0** | 30 | **90% (27/30)** | **0.703** | −0.120 | Silent (26) | **c0: 28/30** |
| 1 | 15 | 6.7% (1/15) | 0.628 | +0.048 | Cushing (14) | c3: 15/15 |
| 2 | 41 | 9.8% (4/41) | 0.534 | −0.035 | Acromegaly/Hyperthyroidism | c4: 41/41 |

→ Neou 单独 86 例样本即可从零复现出三类结构：
1. **NFPT-高甲基化簇**（cluster 0），几乎完全映射到跨队列 NFPT 富集簇 c0（SF1-lineage 类比，非 Møller k4）；
2. **Cushing/TPIT 簇**（cluster 1），映射到 c3；
3. **PIT1 分泌型簇**（cluster 2），映射到 c4。

这与跨队列质心投影结果 1:1 吻合（c0:28, c3:17, c4:41），证明甲基化分型并非由 GSE283928 质心“训练”而来，而是数据驱动的、可在独立 EPIC 队列中自发生成的生物学结构。

---

## 4. 科学意义

1. **甲基化分型框架在第三个独立 EPIC 队列中成立**：Neou 与 GSE283928、GSE207937 共享同一 PCA 甲基化流形，NFPT 与功能性肿瘤的簇映射在两个独立队列中一致。
2. **NFPT 生物学信号跨队列可复制**：NFPT → cluster 0 的投影模式 + Møller 关键基因高甲基化，在 Neou（n=32 NFPT）和 GSE207937（n=45 NS）中同时成立。
3. **跨组学一致性增强**：Møller 关键基因“NFPT 中转录下调 + 甲基化上调”在 Neou 中再次得到支持，强化了这些基因作为 NF-PitNET 甲基化-转录偶联核心锚点的地位。
4. **方法学价值**：整个流程从 idat 到 β、到跨队列投影、到谱系/Møller 评分全部在本地完成，可作为后续外部验证的标准管线。

---

## 5. 诚实边界

- Neou 队列**无结构化术后再生长 / PFS 标签**（clinical history 多为 NA，regrowth 仅零散 free-text），因此本里程碑仍是“甲基化框架捕获 NFPT 生物学”的 surrogate 验证，**非直接预后模型验证**。
- GSE283928 本地矩阵无正常垂体 β，需以论文原文的 118 肿瘤 + 118 正常配对设计作为外部信息，不可在本文件内复现正常对照。
- 谱系代理评分方向与简单预期不一致，已在文中明确声明其局限性。

---

## 6. 产物

- `results/Neou/neou_crosscohort_assignment.csv`
- `results/Neou/gse207937_crosscohort_assignment.csv`
- `results/Neou/neou_lineage_scores.csv`
- `results/Neou/joint_pca_coords.csv`
- `results/Neou/neou_crosscohort_summary.json`
- `results/Neou/neou_crosscohort_figure.png` / `.pdf`
- `results/Neou/neou_internal_clusters.csv`
- `results/Neou/neou_internal_silhouette.json`
- `results/Neou/neou_internal_cluster_summary.json`
- `results/Neou/neou_internal_cluster_figure.png`
- `scripts/neou_export_beta.R`
- `scripts/neou_crosscohort.py`
- `scripts/neou_figure.py`
- `scripts/neou_internal_cluster.py`
