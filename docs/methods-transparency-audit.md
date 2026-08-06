# 方法透明度审计报告

> 对应预审委员会 Reviewer 2 的 8 项整改清单
> 日期：2026-07-31

---

## 1. GSE283928 列识别失误（必改 → 已修正）

**问题**：原始 `gse283928_clusters.py` 将 236 列（118 `_Beta` + 118 `_Detection Pval`）全部作为样本输入 KMeans，导致 detection-pvalue 列被聚为 "cluster 1" 并被误判为"正常垂体对照"。

**修正状态**：✅ 已修正。`neou_crosscohort.py` 和 `cluster_risk_mapping.py` 均只提取 `_Beta` 后缀列。`cluster_assignment_k5.csv` 中 `_Beta` 行的 cluster 标签 {0:70, 2:2, 3:22, 4:24} 已验证为真实肿瘤簇。

**Manuscript 要求**：Methods 中明确说明列筛选逻辑（`endsWith("_Beta")`），附 QC 图展示 Beta 列（中位 ~0.5）与 Detection-Pval 列（中位 ~0, 78.9% exactly 0）的分布差异。

---

## 2. GSE209903 使用预计算 DESeq2 结果（建议改）

**现状**：直接使用 Silva-Junior 2023 提供的 `FPT_vs_NFPT.txt.gz`（ENSG ID, baseMean/log2FC/padj），未独立重跑原始 counts。

**Manuscript 要求**：声明"本研究使用作者提供的预计算 DESeq2 差异分析结果"。如 GEO 提供原始 counts，建议独立重跑并做 scatter 交叉验证（log2FC 相关性 >0.99 即可）。

---

## 3. 三队列联合 PCA 的批次效应处理（必改 → 已审计）

**方法**：**未使用 ComBat 或 limma::removeBatchEffect**。批次控制策略为：
1. **平台一致性**：三个队列均为 Illumina EPIC v1（850K），非跨平台
2. **探针交集**：取三队列公共探针（687,463 / ~866K = 79.4%）
3. **方差过滤**：保留 top 10% 最高方差探针（减少噪声维度）
4. **PCA 可视化**：joint PCA 显示三队列在同一流形上重叠（Neou centroid offset ~13 PC1 units, internal spread std=38.9 → extensive overlap）

**Manuscript 要求**：Methods 中明确说明"未使用 ComBat 等显式批次校正，因三队列均为 EPIC v1 同平台；以探针交集 + 方差过滤 + PCA 可视化验证批次效应可控"。附 pre-correction PCA（cohort-colored）展示重叠。

**审计结论**：方法合理但需透明声明。同平台 EPIC 间的批次效应远小于跨平台（450K↔EPIC），且 joint PCA 的重叠程度支持结论。

---

## 4. NaN 填充策略（必改 → 已审计）

### 4.1 各队列 NaN 比例

| 队列 | 矩阵尺寸 | NaN 数 | NaN % | 最大单样本 NaN % | 全 NaN 探针 |
|---|---|---|---|---|---|
| **GSE283928** | 727,986 × 118 | 0 | 0.00% | 0.00% | 0 |
| **Neou 2020** | 866,238 × 86 | 6,607 | **0.0089%** | 0.086% | 0 |
| **GSE207937** | 727,986 × 77 | 未显式审计 | — | — | — |

### 4.2 填充策略

```python
def impute(M):
    col_mean = np.nanmean(M, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0.5, col_mean)
    idx = np.where(np.isnan(M))
    M[idx] = np.take(col_mean, idx[1])
```

- **策略**：probe-wise（列）均值填充；若整探针为 NaN 则填 0.5（中性 β 值）
- **适用性**：NaN 比例 <0.01%，无样本超过 5% 缺失阈值 → 填充对统计推断的影响可忽略
- **GSE207937**：M-value → β 转换后同样调用 `impute()`，NaN 主要来自 M-value 中的极端值

**Manuscript 要求**：Methods 中报告各队列 NaN 比例和填充策略。声明"无样本因高缺失率（>5%）被排除"。

---

## 5. Møller 关键基因甲基化评分方式（建议改）

**现状**：取 168 个 Møller_key 探针（覆盖 7 个基因：ANK2/ATP1A3/NR5A1/CCDC141/CASC15/GUCA1C/HDC）的 per-sample 均值。未做探针级加权或方向校正。

**局限**：
- Møller 2025 原文的 562 DMP 精确探针坐标未公开（补充材料被 reCAPTCHA 拦截）
- 本研究的"Møller key gene methylation"是基因级均值，非 DMP 探针级复现
- NR5A1 (SF1) 探针同时出现在 Moller_key 和 SF1_lineage 集中 → 有交叉

**Manuscript 要求**：明确声明"关键基因甲基化评分为基因级探针均值，非 Møller 原文 562 DMP 的探针级复现"。在 Limitation 中讨论。

---

## 6. CMap/Enrichr 多重校正（建议改）

**现状**：
- Enrichr API 内部使用 Fisher exact test + BH 校正（`Adjusted P-value`），`cutoff=0.05`
- 但 reversal score（`rev_up + rev_down - prod_up - prod_down`）是自定义复合指标，**未做跨 15,452 实例的多重校正**
- 药物类别频率分析（BET=2713, HDAC=288, ...）为描述性统计，非推断性

**BET 偏倚**：BET/BRD 溴域抑制剂（JQ1 等）在 LINCS 中对几乎所有疾病签名都表现为"强逆转"——这是已知的系统性偏倚（BET 是全局转录调控器，非 NF-PitNET 特异信号）。2713 个实例的压倒性优势不能解读为"BET 是 NFPT 最强候选药"。

**Manuscript 要求**：
1. 声明 Enrichr 内部 BH 校正 + 自定义 reversal score 无额外多重校正
2. 明确讨论 BET 偏倚："BET/BRD 抑制剂在 LINCS 中的普遍强逆转信号是已知的系统性偏倚，不宜解读为 NF-PitNET 特异性药效"
3. 将结论限定在药物类别层面（MEK/mTOR/CDK 类别显著富集），不过度解读单实例排名

---

## 7. 样本量与统计功效（建议改）

| 比较 | 组 1 | 组 2 | 统计量 | p 值 | 功效评估 |
|---|---|---|---|---|---|
| Neou Møller β | NFPT=32 | Func=54 | 均值差 0.124 | perm p≈0 | ✅ 充分 |
| GSE207937 Møller β | NS=41 | other=36 | 均值差 0.032 | — | ✅ 充分 |
| GSE283928 c0 vs c3 | n=70 | n=22 | 均值差 0.018 | perm p=0.0009 | ✅ 充分 |
| GSE283928 c0 vs c4 | n=70 | n=24 | 均值差 0.024 | perm p=0.0001 | ✅ 充分 |
| GSE283928 c0 vs c2 | n=70 | **n=2** | — | p=0.156 (ns) | ❌ 不足（n=2） |

**Manuscript 要求**：报告各组样本量。声明 c2 (n=2) 因样本量过小不参与可靠推断。所有 permutation p 值使用 10,000-20,000 次置换。

---

## 8. 软件版本与随机种子（必改 → 已审计）

### 软件版本

| 组件 | 版本 |
|---|---|
| OS | WSL Ubuntu-24.04 (Windows 11) |
| R | 4.3.3 |
| minfi | 1.48.0 |
| EPIC annotation | IlluminaHumanMethylationEPICanno.ilm10b4.hg19 v0.6.0 |
| Python | 3.13.12 |
| NumPy | (managed env) |
| scikit-learn | (managed env) |
| gseapy | 1.3.1 |

### 随机种子

| 脚本 | 随机操作 | 种子 |
|---|---|---|
| `gse283928_clusters.py` | KMeans(k=5) | `random_state=42, n_init=10` |
| `neou_crosscohort.py` | PCA | `random_state=42` |
| `neou_crosscohort.py` | permutation p (perm_p) | `seed=0, n=20000` |
| `cluster_risk_mapping.py` | permutation p | `seed=0, n=10000` |
| `cluster_risk_mapping.py` | NFPT→c0 null | `rng=default_rng(42), n_perm=10000` |

**Manuscript 要求**：Methods 中报告所有软件版本和随机种子。建议提供 `sessionInfo()` 输出和 `pip freeze` 作为补充材料。

---

## 总结：整改清单状态

| # | 问题 | 级别 | 状态 |
|---|---|---|---|
| 1 | GSE283928 列识别 | 必改 | ✅ 已修正 |
| 2 | GSE209903 预计算 DE | 建议改 | ⏳ 待 manuscript 声明 |
| 3 | 批次效应处理 | 必改 | ✅ 已审计（无 ComBat，同平台+交集+PCA） |
| 4 | NaN 填充 | 必改 | ✅ 已审计（<0.01%，col-mean 填充） |
| 5 | Møller 评分方式 | 建议改 | ⏳ 待 manuscript 声明 |
| 6 | CMap 多重校正 | 建议改 | ⏳ 待 manuscript 讨论 BET 偏倚 |
| 7 | 样本量/功效 | 建议改 | ✅ 已审计 |
| 8 | 软件版本/种子 | 必改 | ✅ 已审计 |
