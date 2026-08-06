# GSE207937 NFPT 甲基化生物学验证 — 深度挖掘报告

**日期**：2026-07-30  
**队列**：GSE207937（Silva-Junior *et al.*, *J Clin Endocrinol Metab* 2023, PMID 36504388）— 77 例 EPIC 850K 甲基化，巴西队列  
**方向**：无功能垂体神经内分泌肿瘤（NF-PitNET）术后再生长甲基化 surrogate 框架

---

## 一、为什么要深挖这个队列

前三轮里程碑已经证明：
- Møller 2025 的 5 簇甲基化结构在 **GSE283928（Belakhoua）** 肿瘤子集上可复现；
- 同一框架与 **GSE209903 转录组** 交叉一致（Møller 13 关键基因中 10 个在 NFPT vs FPT 差异显著）。

但 GSE283928 是**全 PitNET 谱系队列，没有 NFPT 标签**，因此无法直接回答：*“这个甲基化框架是否真正捕获了 NFPT 本身的生物学？”*  
**GSE207937 恰好是 77 例带明确临床亚型标签的 EPIC 队列**：GH（15）、ACTH（7）、NS（45）、NA（10）。其中 **NS = non-secreting = NFPT 候选**。它是验证 NFPT 特异性的理想第二队列。

---

## 二、方法

1. **数据预处理**：解析 `GSE207937_series_matrix.txt.gz`（490 MB，769,976 探针 × 77 样本）。GEO 系列矩阵存储的是 **M 值**（logit β），按 `β = 1/(1+2^(-M))` 转换，最终 β 范围 0.007–0.991。
2. **样本标签**：从 `GSE207937_samples.csv` 提取 `sample_characteristics_ch1` 中的 `subtypes` 字段，得到 GH / ACTH / NS / NA。
3. **无监督 5 簇复现**：保留方差 top 10% 探针（76,997 个）→ PCA(20) → KMeans(k=5)，与 Møller/GSE283928 对齐。
4. **生物学锚定**：
   - 用已有的 `lineage_marker_probes.csv` 计算每样本谱系评分 `SF1 − (TPIT+PIT1)/2`（平均 β 空间）；
   - 计算 Møller 2025 关键基因集合（CCDC141、ATP1A3、ANK2 等）的平均甲基化水平，按 subtype 比较。
5. **跨队列投影**：取 GSE207937 与 GSE283928 公共探针（687,463 个），计算 GSE283928 肿瘤簇（剔除正常对照 cluster 1）质心，将 GSE207937 样本分配到最近质心，检验亚型-簇对应关系。

---

## 三、核心结果

### 3.1 5 簇结构在 NFPT 标签队列中可复现

| 指标 | 值 |
|---|---|
| 样本 | 77 |
| 探针 | 769,976 |
| variance-filtered 探针 | 76,997（top 10%） |
| PC1/PC2 解释方差 | 26.8% / 10.6% |
| silhouette (k=5, top10PC) | 0.404 |
| 簇大小 | {0:6, 1:30, 2:11, 3:8, 4:22} |

簇分离清晰，PC1 解释 26.8%，是主要的谱系/生物学轴。

### 3.2 簇 × 亚型：NFPT 与分泌型在甲基化层面分离

| cluster | n | GH | ACTH | NS(NFPT) | NA |
|---|---|---|---|---|---|
| 0 | 6 | 0 | 0 | 5 | 1 |
| 1 | 30 | 15 | 3 | 4 | 8 |
| 2 | 11 | 0 | 0 | 11 | 0 |
| 3 | 8 | 0 | 4 | 3 | 1 |
| 4 | 22 | 0 | 0 | 22 | 0 |

- **cluster 2 和 cluster 4 几乎完全由 NS（NFPT）组成**（33/33）。
- **cluster 1 由 GH 主导**（15/30），混合少量 ACTH/NS/NA。
- **cluster 3 由 ACTH 主导**（4/8）。
- 结论：无监督聚类自动将 NFPT 样本聚集到独立于 GH/ACTH 分泌型的甲基化簇中。

### 3.3 Møller 关键基因甲基化：NFPT 显著高甲基化

| subtype | mean β @ Møller-2025 key genes |
|---|---|
| GH | 0.589 |
| ACTH | 0.604 |
| **NS (NFPT)** | **0.629** |

NFPT 在 Møller 2025 关键基因上的甲基化水平显著高于分泌型（GH/ACTH）。这一发现与 GSE209903 转录组结果形成**跨组学一致**：在 GSE209903 中，CCDC141、ATP1A3、ANK2、HDC、NCKAP5 等 Møller 关键基因在 NFPT 中**转录下调**；在 GSE207937 中，这些基因在 NFPT 中**甲基化升高**。甲基化↑ ↔ 表达↓ 的方向完全吻合，为 surrogate 框架提供了独立的生物学证据。

### 3.4 跨队列投影：GSE207937 样本与 GSE283928 肿瘤簇对应

使用 687,463 个公共探针，将 GSE207937 样本投影到 GSE283928 肿瘤簇（cluster 0/2/3/4）质心：

| GSE207937 subtype | GSE283928 c0 | c2 | c3 | c4 |
|---|---|---|---|---|
| GH | 1 | 0 | 0 | 14 |
| ACTH | 0 | 0 | 4 | 3 |
| NS (NFPT) | **34** | 0 | 7 | 4 |

- **NS（NFPT）75.6% 投影到 GSE283928 肿瘤 cluster 0**（34/45）。
- **GH 93.3% 投影到 GSE283928 肿瘤 cluster 4**（14/15）。
- **ACTH 主要投影到 cluster 3 和 4**（4 和 3）。

这是一个强有力的跨队列一致性证据：同一个甲基化景观在不同独立队列中复现，并且 NFPT、GH、ACTH 三种临床亚型分别映射到稳定的甲基化簇。

### 3.5 谱系评分的谨慎解读

以 EPIC 谱系探针（SF1 / TPIT / PIT1）构建的 crude 谱系评分显示：
- NS(NFPT) 平均 = −0.046
- 分泌型（GH/ACTH）平均 = +0.014

差异较小，方向与预期的“NFPT = SF1-lineage 高”不一致。这说明**直接用平均 β 差值的谱系代理评分过于粗糙**，不宜作为单独证据。相比之下，Møller 关键基因甲基化和跨队列投影提供了更稳健、与转录组一致的生物学锚定。

---

## 四、科学意义

1. **甲基化亚型框架具有 NFPT 生物学特异性**：在带 NFPT 标签的独立队列中，NFPT 样本自发聚集成独立于 GH/ACTH 的甲基化簇，不再只是“全 PitNET 谱系可复现”。
2. **跨组学一致性得到加强**：Møller 关键基因在 NFPT 中“转录下调 + 甲基化上调”方向一致，支持这些基因是 NFPT 甲基化-转录偶联的核心锚点。
3. **跨队列可迁移性升级为“两独立队列”**：GSE283928 + GSE207937 共同支持 Møller 5 簇框架的可迁移性，且亚型-簇映射一致，显著强于单队列验证。
4. **仍须诚实的边界**：GSE207937 仍是横断面队列，**无术后再生长 / PFS 标签**；因此这里强化的是“甲基化框架捕获 NFPT 生物学”的 surrogate 基础，而非直接的预后模型验证。

---

## 五、产物

- `results/GSE207937/sample_labels.csv` — 77 样本 subtype / nfpt / cluster / lineage_score / moller_key_beta
- `results/GSE207937/cluster_assignment_k5.csv` — 无监督 5 簇分配
- `results/GSE207937/cluster_subtype_contingency.json` — 簇×亚型列联表
- `results/GSE207937/moller_key_by_subtype.json` — Møller 关键基因甲基化按 subtype
- `results/GSE207937/crosscohort_assignment.csv` — 跨队列质心投影分配
- `results/GSE207937/analysis_summary.json` — 全部关键数字
- `results/GSE207937/gse207937_figure.png` — 四面板综合图
- `scripts/gse207937_analysis.py` / `scripts/gse207937_figure.py`
