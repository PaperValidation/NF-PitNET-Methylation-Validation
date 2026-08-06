# 最终稿全面检查报告（含 Table 与 Figure）

**检查对象**: `docs/manuscript-v1.md`（ARS 第一版初稿）+ 全部配套表格/图/补充材料
**检查日期**: 2026-08-04
**检查方式**: ARS academic-paper-reviewer 工作流 → 逐数字对照结果文件（JSON/CSV/PNG 实体）

---

## 总结论

稿件科学内容与结果文件**高度一致**（18 个关键数字组全部核对通过，引用无孤儿），但存在 **2 个 CRITICAL**（表格实体缺失、图编号跳号）、**4 个 MAJOR**（数字口径、图合并、元数据失实）、**8 个 MINOR**。投稿前必须解决 CRITICAL。

---

## 一、CRITICAL（必须修复）

### C1. Table 1 和 Table 2 只有标题，无实体内容
- **位置**: `manuscript-v1.md` L307-309（Table 1 标题下无表格；Table 2 标题下无表格）
- **正文引用**: §2.1 引用 Table 1（"summarized in Table 1"）、§3.2 引用 Table 2（"**Figure 1B**; **Table 2**"）
- **问题**: 正文引用了不存在的表格。Table 1（cohort characteristics）和 Table 2（cross-cohort projection contingency）均无实际数据表。
- **修复建议**: 从 `results/Neou/neou_crosscohort_summary.json`（stage/lineage counts、nfpt_vs_cluster、threeway_cluster_counts）、`results/GSE207937/cluster_subtype_contingency.json`、`results/Moller2025/strict_projection/moller_strict_assignment.csv`、`results/Moller2025/gse115783_validation/gse115783_assignment.csv` 生成两张真实表格。

### C2. 主图编号跳号：无 Figure 3
- **位置**: `manuscript-v1.md` §7 图注（Figure 1 → Figure 2 → **Figure 4** → Figure 5 → Figure 6）；正文同样引用 Figure 4/5/6
- **原因**: 原 Figure 3（CMap 图）在 nature-reviewer 整改时移至 Supplementary Figure S3，但主图编号未重排
- **问题**: 投稿稿图必须连续编号（Figure 1-5）。当前从 2 直接跳 4。
- **修复建议**: 重排编号 Figure 4→3、Figure 5→4、Figure 6→5，并同步正文所有引用（§3.3/§3.5/§3.8/§3.9/§4/§5 等）；或保留 Figure 3 位置并补充说明。推荐重排。

---

## 二、MAJOR（建议修复）

### M1. NN-mixing 数字存在两套口径，稿件未说明
- **位置**: §3.1（L153）、§2.3（L114）
- **稿件表述**: "38.9% of nearest neighbors... Neou observed 0.51 (74% of its 0.69 random expectation), GSE207937 observed 0.51 (70% of its 0.73), GSE283928 observed 0.22 vs 0.58 expected... observed 0.389... 置换零 0.66"
- **实测对照**:
  - `results/Neou/batch_effect_neighbor_fraction.json`: overall 0.3893, per-cohort {GSE207937 0.5065, GSE283928 0.2243, Neou 0.5109}, random_per_cohort {GSE207937 0.726, GSE283928 0.580, Neou 0.694} → **稿件数字与这个文件完全吻合**（0.51/0.51/0.22 ✓）
  - `results/Neou/nn_mixing_permutation.json`（新脚本）: observed_overall_mixing **0.33**, per_cohort {GSE207937 0.4476, GSE283928 0.2051, Neou 0.3961}, size_weighted_random_expectation 0.6549, perm_null_mean 0.6572, p=2e-4
- **问题**: 同一概念出现两个 observed（0.389 vs 0.33）和两套 per-cohort 数字。差异源于两脚本 k 邻近处理细节（是否排除自身/距离度量）。摘要记忆已注明"原稿 0.389 保留"，但稿件正文未透明披露双口径。
- **修复建议**: 在 §2.3 方法中注明"0.389 为主分析口径（batch_effect_neighbor_fraction.py），0.33 为独立复算（nn_mixing_permutation.py，k=15 排除自身），两者同量级、置换 p<0.001 结论一致"；或统一到单一脚本重跑后全文更新。

### M2. Figure 5 实际是三个独立 PNG，未合并为单图
- **位置**: §7 图注 Figure 5（描述为 5A/5B/5C 单图）
- **实际文件**: `results/Moller2025/strict_projection/km_strict_riskhi.png`（900×675）+ `forest_strict.png`（900×570）+ `fig5c_lineage_bar.png`（825×570）
- **问题**: 投稿需单张多面板图（Figure 5 含 A/B/C 三面板）。当前三图分辨率偏低（~900px 宽，300dpi 需 ~2500px）。
- **修复建议**: 用 matplotlib 合并三面板为单图并提高 dpi（≥300），重导出。

### M3. Draft Metadata 与实际文件不符
- **位置**: §Draft Metadata（L374-392）
- **不符项**:
  - "Total Word Count ~5,700" → 实测 **5,887**（+17% 不是 +14%）
  - "Em dash density: Reduced from 64 to ~20" → 实测 v1 **em dash = 0**（全部已替换）
  - "Abstract ~177" → 摘要实测约 177 词 ✓（此项 OK）
- **修复建议**: 更新为实测值 5,887 / 0 em dash / 一致性声明。

### M4. Figure 6 图注颜色标注无法从文件确认 + Figure 分辨率不足
- **位置**: §7 Figure 6 图注（"GSE115783 NFPA tumors, orange"）
- **实际文件**: `results/Moller2025/external_validation_consolidated/external_validation_bar.png`（1425×780）
- **问题**: 无法读图确认颜色分组是否与图注一致（需人工目检）；1425px 宽低于 300dpi 投稿要求（183mm 宽约需 2160px）。
- **修复建议**: 人工目检 external_validation_bar.png 是否区分"EPIC 三队列 vs HM450"颜色；提高导出 dpi。

---

## 三、MINOR（建议微调）

| # | 位置 | 问题 | 建议 |
|---|------|------|------|
| m1 | 摘要（L52） | "is often misread as 'aggressive NF-PitNET'" 无引用断言，与 §4 已改为 "tempting to equate...is incorrect" 的措辞不一致 | 摘要同步改为 "is sometimes conflated with 'aggressive NF-PitNET'"（弱化断言） |
| m2 | Conclusion（L261） | "A NF-PitNET-associated" 语法错误（NF 读元音开头） | 改为 "An NF-PitNET-associated" |
| m3 | §7 Figure 5C 图注 | "101/117 (86%)" vs 正文 "101/117 (86.3%)" | 统一为 86.3% |
| m4 | §2.4 vs §3.1 | 方法写 20,000 permutations（assignment null），结果写 5,000 shuffles（mixing null），未区分两套置换 | §2.4 注明两种 null 及其置换次数 |
| m5 | Data Availability | "deposited in scripts/ and results/ directories of the project repository" 无实际 repo URL | 投稿前补真实公开仓库链接 |
| m6 | §Supplementary Materials | 指向本地路径（`docs/supplementary-*.md`、`results/...`） | 投稿时替换为 Figshare/Zenodo 或期刊补充文件编号 |
| m7 | §3.9 | "Radiotherapy was excluded because it showed complete separation" — rt_yn 在 JSON 中为 ConvergenceError，方法需注明排除理由 | 补充一句"RT 与事件完全分离导致 Cox 不收敛" |
| m8 | §2.6/§3.4 | "NR5A1/SF1 promoter hypermethylation is a generic PitNET-wide mark (ref 7)" — 正文两处重复表述 | 保留一处，另一处简化为交叉引用 |

---

## 四、核对通过项（数字与结果文件一致）

以下 18 组关键数字**全部核对通过**，无错误：

| 类别 | 数字 | 状态 |
|------|------|------|
| 四队列 c0 映射 | 81.2/75.6/86.3/94.1%（均值 84.3%，EPIC-only 81.0%） | ✅ |
| GSE115783 | 32/34=94.1%，正常对照 2/6=33.3% | ✅ |
| 预后 HR | k3/k4/k5 vs k1/k2: HR 2.14 (0.94-4.89, p=0.072)，logrank 0.065 | ✅ |
| 校正 HR | 2.09 (0.90-4.86, p=0.088)；术前体积 1.05/cc (p=5.3e-3) | ✅ |
| k4-vs-rest | HR 1.82 (0.76-4.37, p=0.18) | ✅ |
| 事件率 | 28/111=25.2%；k1 22.2/k2 13.3/k3 30.8/k5 27.8/k4 36.8% | ✅ |
| 562-DMP 方向一致性 | 0.621 (p=8.2e-8) / 0.566 (p=1.1e-3) / 0.525 (p=0.13) | ✅ |
| k4→c3/c4 | 1.000, p=1.6e-30；GSE207937 secreting 1.000, p=7.7e-34；Neou int1 0.591, p=0.035 | ✅ |
| SF1 簇映射 | k2→c0 0.982 (p=1.2e-30)、k5→c0 0.830、k3→c0 0.709、k1→c0 0.537 (ns) | ✅ |
| 562-DMP LOO | 99.1% (110/111)、RF 93.7% (104/111, κ=0.92)、majority 27.0% | ✅ |
| 归一化 | 0.326/0.326/0.337/0.314，一致 100/100/98.8/98.8% | ✅ |
| PCA | PC1 29.2%、PC2 9.7%、cum20 63.4% | ✅ |
| 内部聚类 | k=3 (silhouette 0.4045)；int0 n=30/90% NFPT；28/30→c0、15/15→c3、41/41→c4 | ✅ |
| strict 117 | 101/117=86.3%；c0 dist 1153 < c3 1235 < c4 1400；c0-vs-non HR 1.63 (p=0.30) | ✅ |
| 探针数 | 687,463 common；68,746 top-var；86,623 Neou-only；333,574 HM450；725,762 Moller | ✅ |
| CMap | n=15,452 reversal-positive；BET/BRD 2,713；Top-20 表与 drug_reversal_top40.csv 一致 | ✅ |
| 关键基因 | 168 CpG；ANK2/ATP1A3/NR5A1/CCDC141/CASC15/GUCA1C 注释正确 | ✅ |
| 引用 | 1-18 全部正文引用且全部列出，**无孤儿引用、无重复编号** | ✅ |

---

## 五、Table/Figure 实体清点

| 编号 | 图注标题 | 实体文件 | 状态 |
|------|---------|---------|------|
| Figure 1 | Cross-cohort validation | `results/Neou/neou_crosscohort_figure.png` (1950×1650) | ✅ 存在 |
| Figure 2 | Internal clustering | `results/Neou/neou_internal_cluster_figure.png` (1650×1350) | ✅ 存在 |
| ~~Figure 3~~ | （CMap，已移 S3） | — | ⚠️ 编号跳号（C2） |
| Figure 4 | 562-DMP probe validation | `results/Moller2025/moller_dmp_validation_figure.png` (3197×1898) | ✅ 存在 |
| Figure 5 | KM + forest + lineage bar | 三个独立文件（M2） | ⚠️ 需合并 |
| Figure 6 | External validation bar | `results/Moller2025/external_validation_consolidated/external_validation_bar.png` (1425×780) | ✅ 存在（分辨率偏低） |
| Fig S3 | Drug candidates | `results/CMap/supplementary_figure_s3.png` (2089×1076) | ✅ 存在 |
| Table 1 | Cohort characteristics | **无实体** | ❌ C1 |
| Table 2 | Projection contingency | **无实体** | ❌ C1 |
| Table 3 | Normalization sensitivity | v1 L313-318 内联表 ✅ | ✅ 与 ST3 CSV 一致 |
| ST1/ST2/ST3 | 补充表 | `results/Neou/supplementary_table*.csv` | ✅ 存在 |
| Supp Methods | S1+S2 | `docs/supplementary-methods.md` | ✅ 存在 |

---

## 六、修复优先级建议

1. **[C1] 生成 Table 1/2 实体**（数据齐全，脚本 30 分钟可完成）
2. **[C2] 重排图编号** Figure 4→3、5→4、6→5 + 正文引用同步
3. **[M2] 合并 Figure 5 三面板 + 提高 dpi**
4. **[M1] 双口径 NN-mixing 透明化**（方法中注明两脚本差异）
5. **[M3] 更新 Draft Metadata**（5,887 词 / 0 em dash）
6. MINOR 项逐条微调（m1-m8，其中 m2 语法为必改）

---

*报告生成：ARS academic-paper-reviewer 工作流，基于结果文件逐项核对，非记忆推断。*
