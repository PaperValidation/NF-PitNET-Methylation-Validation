# 写作团队交接文件 — NF-PitNET 甲基化分类学跨队列验证

> **交接方**：Biomni 生物医学专家团（分析侧）
> **接收方**：写作类专家团（如 `academic-research-team`）
> **生成日期**：2026-10-18
> **分析状态**：✅ 全部完成并经数字级复核，逐项可追溯
> **交接范围**：分析产物、证据链、图表素材、诚实边界。**本团不承接任何论文成文工作。**

---

## 0. 一句话结论（供写作定位）

非功能型垂体神经内分泌肿瘤（NF-PitNET）在**四个独立队列、两个芯片平台**上稳定映射到 GSE283928 的 **c0 簇（SF1-lineage 甲基化身份）**；而文献中被视为"最高风险"的 **k4（TPIT/PIT1, HR 3.63）并不对应 c0**，而是映射到分泌型簇 c3/c4（探针级 p=1.6×10⁻³⁰）。因此"NF-PitNET = 最高风险 k4"是一个流行误读。真正承载再生长风险的是 c0 **内部**的谱系内细分（k3/k4/k5 轴）。

---

## 1. 核心资产位置

| 资产 | 路径 | 说明 |
|---|---|---|
| **终稿** | `docs/manuscript-v2.md` | 已按 *Neuro-Oncology* 要求校准，68KB，16 条参考文献 |
| Cover letter 草稿 | `docs/cover-letter.md` | 待填作者信息 |
| 投稿前终检清单 | `docs/submission-checklist.md` | A–F 六组，可逐项勾选 |
| 待填信息清单 | `docs/author-action-checklist.md` | 6 处 `[AUTHOR_ACTION]` 逐条列明 |
| 数据来源权威说明 | `docs/data-availability.md` | 6 队列来源/DOI/本地布局 |
| 分析产物 | `results/`（75 文件，17 图件） | 见第 3、4 节 |
| 分析脚本 | `scripts/`（54 py + 18 R） | 全部经 `py_compile` 通过 |

---

## 2. 数字 → 证据来源 完整对照表（写作时每个断言都必须能落到这里）

> **使用规则**：正文中出现的任何数字，必须能在下表找到"产物路径"。若写作过程中需要新数字，请勿自行推算，回到本团重新出数。

### 2.1 四队列外部验证（Figure 5 / Table 2 / 摘要）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| Neou 2020 NFPT → c0 | 26/32 = **81.2%** | `results/Moller2025/external_validation_consolidated/external_validation_consolidated.json` |
| GSE207937 NS → c0 | 34/45 = **75.6%** | 同上 |
| Møller 原始 idat → c0 | 101/117 = **86.3%** | 同上 + `results/Moller2025/strict_projection/strict_projection_summary.json` |
| GSE115783 NFPA → c0 | 32/34 = **94.1%** | 同上 + `results/Moller2025/gse115783_validation/gse115783_summary.json` |
| 四队列均值 | **84.3%**；EPIC-only 均值 **81.0%** | 同上（json 内含 mean 字段） |
| GSE115783 正常垂体对照 → c0 | 2/6 = **33.3%**（4/6 → c4） | `gse115783_summary.json` |
| Møller 投射平均距离 | c0 **1153** < c3 1235 < c4 1400（c2 1568） | `strict_projection_summary.json` → `centroid_mean_distance` |
| Møller 投射共有 cg 探针 | **725,762** | 稿件 §3.8 明载（由 `moller_strict_projection.py` 计算） |
| GSE115783 跨平台桥接探针 | **333,574** | `gse115783_summary.json` → `common_probes` |

### 2.2 决定性证据：k4 ≠ c0（Figure 3 / 摘要 / 讨论）
> ⚠️ **口径警告（最容易搞错）**：稿件引用的是 `moller_k_cluster_mapping.json`（**k → 具体 cluster** 口径）。另有一个 `moller_dmp_validation_summary.json`（**high-risk vs low-risk 组总体**口径，其中 k4 frac=0.0, p=1.0）是**不同分析**，**不可**用它否定下表结论。两个文件同时存在于结果目录，请勿混引。

| 断言 | 数值 | 产物路径 |
|---|---|---|
| **k4 → GSE283928 c3/c4** | frac **1.000**, p=**1.6e-30**（两簇各自 99/99） | `results/Moller2025/moller_k_cluster_mapping.json` |
| k4 → c0 | frac **0.0**（99 个探针全反向） | 同上 |
| k4 → GSE207937 secreting | **1.000**, p=**7.7e-34** | 同上 |
| k4 → GSE207937 NS | **0.0** | 同上 |
| k4 → Neou int1 (Cushing/TPIT) | 0.591, p=**0.035** | 同上 |
| k2 → c0 | 0.982, p=**1.2e-30** | 同上 |
| k5 → c0 | 0.830, p=**6.5e-12** | 同上 |
| k3 → c0 | 0.709, p=**1.3e-4** | 同上 |
| k1 → c0 | 0.537, **ns** (p=0.29) | 同上 |
| k2 → GSE207937 NS | 0.813, p=**6.1e-13** | 同上 |
| k5 → GSE207937 NS | 0.778, p=**6.1e-10** | 同上 |
| k1 → GSE207937 NS | 0.654, p=**9.2e-4** | 同上 |
| 562-DMP 总数 / 各簇归属 | 562；k1=108, k2=124, k3=102, k4=110, k5=118（union 562，三队列共有 470） | 同上 + `moller_dmp_validation_summary.json` |

### 2.3 谱系内细分可复现性（Figure 4 / §3.8）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| 562-DMP 最近质心 LOO | **99.1%**（110/111） | `results/Moller2025/within_lineage_562dmp/moller_562dmp_summary.json` |
| 562-DMP 随机森林 LOO | **93.7%**（104/111），**κ=0.92** | 同上 |
| 多数类基线 | **27.0%**（k2 为多数类） | 同上 |
| 错分位置 | 全部发生在相邻 SF1 簇 k2↔k3 之间；k4 从未越界 | 同上 `nearest_centroid_LOO_confusion` |
| 样本数 / k 分布 | n=111；k1=18, k2=30, k3=26, k4=19, k5=18 | 同上 |

### 2.4 直接预后验证（Figure 4A/4B）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| 队列 / 事件 | n=**111**（117 例中 6 例无随访排除）；事件 28；事件率 25.2% | `results/Moller2025/prognostic_validation/prognostic_validation_summary.json` |
| **k3/k4/k5 vs k1/k2** Cox | HR **2.14**, 95%CI **0.94–4.89**, p=**0.0715** | 同上 `high_vs_low` |
| log-rank | p=**0.0652** | 同上 |
| 调整后（年龄/性别/侵袭/术前体积） | HR **2.09**, 95%CI **0.90–4.86**, p=**0.0884** | `strict_projection_summary.json` → `multivariable.risk_hi`；`prognostic_validation_multivariable.json` |
| 累积再生长 Kruskal-Wallis | p=**8.5e-3**（0.00849） | `prognostic_validation_summary.json` → `regrow_pct_kruskal_p` |
| 术前体积（独立显著） | HR **1.05–1.07/cc**, p=**5.3e-3**（单变量 p=1.2e-4, HR 1.069） | `prognostic_validation_multivariable.json` → `univariable.preop_vol` |
| 侵袭性生长（单变量） | HR 2.43, p=0.034 | 同上 `univariable.invasive_bin` |
| k4 vs rest（单独） | HR 1.82, p=**0.181**（不显著） | `prognostic_validation_summary.json` → `k4_vs_rest` |
| 各簇事件率 | k1 22.2% / k2 13.3% / k3 30.8% / k4 36.8% / k5 27.8% | 同上 `per_k` |
| 各簇中位再生长% | k1 21.2 / k2 20.6 / k3 59.1 / k4 55.9 / k5 74.1 | 同上 |
| **c0 vs non-c0 预后（阴性结果）** | HR **1.63**, p=**0.30**；log-rank p=0.304 | `strict_projection_summary.json` → `strict_nonc0_vs_c0` |

> **写作要点**：c0 粗分**不能**预测预后这一阴性结果必须如实写出——它恰恰证明"风险在谱系内细分而非谱系身份本身"，是本文逻辑闭环的关键一环，不是失败。

### 2.5 Neou 内部从头聚类（Figure 2）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| 最优 k | **k=3** | `results/Neou/neou_internal_cluster_summary.json` → `best_k` |
| silhouette by k | k2 0.4026 / **k3 0.4045** / k4 0.3256 / k5 0.3023 / k6 0.2921 | `results/Neou/supplementary_tables_summary.json` → `st2_silhouette_by_k` |
| int0 | n=30，NFPT 27（**90.0%**），dominant SF1，Møller-key β 0.7027，→ c0 (28) / c3 (2) | `supplementary_tables_summary.json` → `st2` |
| int1 | n=15，NFPT 1（6.7%），dominant TPIT，CushingOvert，→ c3 (15/15) | 同上 |
| int2 | n=41，NFPT 4（9.8%），dominant PIT1，Acromegaly，→ c4 (41/41) | 同上 |

### 2.6 联合 PCA 与批次效应（Figure 1A）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| 共有探针 | **687,463**（3 EPIC 队列） | `results/Neou/neou_crosscohort_summary.json` |
| 方差过滤后 | 68,746（top 10%） | 同上 |
| PC1 / PC2 解释率 | **29.23% / 9.67%** | 同上 |
| PC1–20 累积 | **63.44%** | 同上 |
| 联合 PCA 样本数 | 281（Neou 86 + GSE283928 118 + GSE207937 77） | `results/Neou/nn_mixing_permutation.json` |
| NN-mixing 观测值 vs 零分布 | 0.330 vs 0.657（null 95%CI 0.638–0.673），**p=2.0e-4** | 同上（perm_n=5000） |
| 跨队列近邻比例 | 总体 0.389；Neou 0.511 / GSE207937 0.506 / GSE283928 0.224 | `results/Neou/batch_effect_neighbor_fraction.json` |

> 这两项是"三队列共享同一甲基化流形、非批次驱动"的量化证据。

### 2.7 归一化敏感性（Table 3 / §3.1）
| 归一化 | c0/86 | c0 占比 | 与 preprocessRaw 一致性 | 
|---|---|---|---|
| preprocessRaw | 28 | 0.326 | 基准 |
| noob | 28 | 0.326 | **100.0%** |
| funnorm | 29 | 0.337 | **98.8%** |
| quantile | 27 | 0.314 | **98.8%** |

产物：`results/Neou/supplementary_tables_summary.json`（`st3`）；重分配样本仅 **P009、P082** 两例；c0 占比区间 **0.314–0.337**（波动 <2.3 个百分点）。

> 注：归一化敏感性分析**未在产物中单独记录探针数**（脚本按 `intersect(rownames(beta_raw), cen$probe)` 动态取交集）。若写作时需要该数字，须回分析侧重跑取数，**不要**把 Møller 投射的 725,762 挪用到此处。

### 2.8 Møller-key 基因甲基化（Figure 1C/1D）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| Neou NFPT vs Func | 0.6873 vs 0.5627（perm p≈0） | `results/Neou/neou_crosscohort_summary.json` |
| GSE207937 NS vs other | 0.6287 vs（GH 0.5889 / ACTH 0.6044） | `results/GSE207937/analysis_summary.json` |
| 谱系代理评分（局限） | Neou NFPT −0.1173 vs Func −0.0104；GSE207937 NS −0.0458 vs secretory +0.0135 | 同上两文件 |

> ⚠️ **必须先读 §5.2**：这个信号主要反映 **SF1-lineage 谱系身份**（探针集含 *NR5A1*/SF1 本身），不能当作 k4 风险信号；crude SF1-score 也不宜作独立证据（NFPT 反而更低）。写作时切勿把它包装成"风险信号"。

### 2.9 药物重定位（Supplementary Figure S3，探索性）
| 断言 | 数值 | 产物路径 |
|---|---|---|
| 分析化合物数 | 15,452 | `results/CMap/cmap_summary.json` |
| 类别计数 | BET/BRD 2713、HSP90 359、HDAC 288、SFK/ABL 239、mTOR/PI3K 223、EGFR/HER 118、CDK 101、MEK 61、RAF/MEK 32、GC 18、Hormone 9、JAK/STAT 8 | 同上 |
| Top 候选 | selumetinib、PD-0325901、dasatinib、WZ-4-145 等 | `results/CMap/drug_reversal_top40.csv` |

---

## 3. 图件清单（全部已生成，`results/` 下）

| 图 | 文件 | 大小 | 对应稿件 |
|---|---|---|---|
| Figure 1 | `results/Neou/neou_crosscohort_figure.png`（+ `.pdf` 矢量版） | 222K | 联合 PCA、跨队列投射、Møller-key |
| Figure 2 | `results/Neou/neou_internal_cluster_figure.png` | 167K | Neou 内部 k=3 聚类 |
| Figure 3 | `results/Moller2025/moller_dmp_validation_figure.png`（+ `.pdf`） | 372K | 562-DMP 方向一致性热图 |
| Figure 4 | `results/Moller2025/figure4_merged.png` | 557K | A: KM / B: forest / C: 独立谱系投射 |
| Figure 4A 单独 | `results/Moller2025/strict_projection/km_strict_riskhi.png` | 46K | KM 曲线 |
| Figure 4B 单独 | `results/Moller2025/strict_projection/forest_strict.png` | 33K | Forest plot |
| Figure 4C 单独 | `results/Moller2025/strict_projection/fig5c_lineage_bar.png` | 32K | 谱系投射柱状 |
| Figure 5 | `results/Moller2025/external_validation_consolidated/external_validation_bar.png` | — | 四队列 c0 映射柱状 |
| Supp Fig S3 | `results/CMap/supplementary_figure_s3.png` | — | 药物重定位类别 |
| 其他备用 | `results/Moller2025/prognostic_validation/{km_moller_taxonomy_regrowth,forest_moller_cox}.png`、`results/Moller2025/within_lineage_562dmp/moller_562dmp_confusion.png`、`results/Moller2025/gse115783_validation/gse115783_lineage_bar.png`、`results/GSE207937/gse207937_figure.png`、`results/CMap/cmap_figure.png` | — | 可作补充图 |

**共 17 个 PNG/PDF 图件。**

> ⚠️ 投稿前需按 *Neuro-Oncology* 规范重新导出（RGB、≥300 dpi 彩色 / ≥600 混合 / ≥1000 线条、字≥7pt Arial、宽≤183mm）。现有图为分析用分辨率，**不能直接投稿**。

---

## 4. 表格（Table 1–3 已在稿件正文；补充表有 CSV 原始数据）

| 表 | 内容 | 原始数据 |
|---|---|---|
| Table 1 | 队列特征（5 队列平台/预处理/N/组成） | 稿件内 |
| Table 2 | 跨队列投射列联表（4 组） | 见 §2.1 产物 |
| Table 3 | 归一化敏感性 | `results/Neou/supplementary_table3_normalization.csv` |
| Supp Table 1 | 联合 PCA 队列质心与方差解释 | `results/Neou/supplementary_table1_pca.csv`、`supplementary_table1_pca_cluster_centroids.csv` |
| Supp Table 2 | Neou 内部聚类特征 | `results/Neou/supplementary_table2_internal_clusters.csv` |
| Supp Table 3 | 归一化敏感性四簇细分 | `results/Neou/supplementary_table3_normalization.csv` |
| 逐样本分配明细 | 四队列归属 | `results/Moller2025/strict_projection/moller_strict_assignment.csv`、`results/Moller2025/within_lineage_562dmp/moller_562dmp_assignment.csv`、`results/Moller2025/gse115783_validation/gse115783_assignment.csv`、`results/Neou/neou_crosscohort_assignment.csv`、`results/Neou/neou_internal_clusters.csv`、`results/GSE207937/crosscohort_assignment.csv` |
| 患者级生存数据 | Møller 111 例 time-event | `results/Moller2025/prognostic_validation/moller_patient_timeevent{,_enriched}.csv` |

补充方法与表格的完整文字版：`docs/supplementary-methods.md`、`docs/supplementary-tables.md`。

---

## 5. 诚实边界（**写作时不可越界**，这些是分析侧明确划定的红线）

### 5.1 预后结论的强度
- HR 2.14 / p=0.072、调整后 HR 2.09 / p=0.088 —— **均为 borderline，不显著**。
- 必须使用"**suggestive**"、"**trend toward**"、"**borderline significance**"等措辞。
- **禁止**写成 "significantly predicts regrowth"、"prognostic biomarker"、"validated prognostic classifier"。
- 唯一显著的是**术前肿瘤体积**（HR 1.05–1.07/cc, p≈5e-3）——这是临床变量，不是甲基化标志物。
- 终点是**再手术/再干预（surrogate）**与影像再生长，**不是 PFS/OS**；公共队列**无样本级 PFS 标签**。
- 样本量 n=111、事件 28 —— 功效有限，必须在 Limitations 里点明。

### 5.2 SF1-lineage 作为药敏标志物 = 假说，零直接证据
- **"谱系身份 → 药敏"目前没有任何直接证据**。药物重定位是**探索性、假设生成性**（hypothesis-generating），仅基于转录签名反演。
- 候选药证据层级（务必按此分级表述，不可抬高）：
  - **Tier 1（有真实 PitNET 证据）**：HDAC/panobinostat（Acta Neuropathol Commun 2024）≈ mTOR/everolimus（Pituitary 2023）> BET/BRD4（Neuro-Oncol 2020）
  - **Tier 2（低置信，需功能验证）**：MEK（仅促肾上腺皮质激素型有证据）、CDK4/6（仅 1 例 NFPA 个案）
- **禁止**暗示 MEK/CDK4/6 已有临床依据；**禁止**把 SF1-lineage 写成药敏预测标志物。
- 文献校准详版：`docs/deep-research-nf-pitnet-drug-lineage-validation.md`；原代细胞验证路线图：`docs/nfpa-primary-cell-functional-validation-plan.md`（尚未做湿实验）。

### 5.3 曾出现过的错误结论（务必不要复活）
- ❌ "c0 = Møller 最高风险 k4" —— 已被 562-DMP 探针级证据**证伪**。
- ❌ "GSE283928 cluster 1 = 正常垂体对照" —— 那是 detection-pvalue 伪影列。
- ❌ 把 Møller-key 基因高甲基化当作"k4 风险信号" —— 它反映 SF1 谱系身份。
- 这些订正过程见 `docs/moller-dmp-probe-validation.md`、`docs/narrative-reframing-options.md`。

### 5.4 c0 粗分的阴性结果必须保留
c0 vs non-c0 预后 HR 1.63 / p=0.30 —— **必须如实报告**。它证明"风险轴是谱系内细分"，是支撑主论点的关键，删掉会破坏逻辑闭环。

---

## 6. 参考文献状态

- **16 条，Vancouver 编号制，1–16 连续，全部经 PubMed/Crossref 核验为真实存在，且 16 条全部被正文引用**（无孤儿引用、无未被引用的条目）。
- 已删除 2 条**编造文献**（原 [2] Ntali/Grossman 2014、原 [4] Sievers 2020）——详见 `docs/reference-verification-report.md`。
- ⚠️ **写作团队若要新增参考文献，必须逐条到 PubMed/Crossref 核实存在性**（作者+期刊+年检索，不能只按题名，题名本身可能是幻觉）。历史上的编造引用正是这样产生的。

---

## 7. 格式合规状态（已程序化校验，改动后需复验）

| 项目 | 要求 | 当前 | 状态 |
|---|---|---|---|
| Structured abstract | ≤250 词 | 202 词 | ✅ |
| Keywords | ≤5 | 5 | ✅ |
| Key Points | 2–3 条，每条 ≤85 字符 | 3 条（79 / 71 / 77） | ✅ |
| Importance of the Study | ≤150 词 | 114 词 | ✅ |
| Running title | ≤50 字符 | 45 字符 | ✅ |
| References | Vancouver 编号 | 16 条，连续 | ✅ |
| AI Disclosure | 期刊强制 | 存在 | ✅ **不可删** |
| Data Availability / Ethics / CRediT / COI / Funding | 均强制 | 前二已写，后三待填 | ⚠️ |

---

## 8. 待用户/写作团队处理（6 处 `[AUTHOR_ACTION]`）

| 位置 | 内容 |
|---|---|
| manuscript-v2.md L18 | 作者列表 + 单位 + ORCID |
| L22 | 通讯作者 + 联系方式 |
| L360–364 | CRediT 贡献（骨架已给） |
| L368 | COI 确认（全体作者） |
| L372 | Funding 声明（无资助则用期刊模板句） |
| Zenodo `10.5281/zenodo.21840024` | creators 占位 "Tidus to be replaced" → 真实作者名单 |

逐条清单见 `docs/author-action-checklist.md`。

---

## 9. 可复现性与环境（**2026-10-18 有变更，务必注意**）

- **磁盘已重组**：原 G: 盘的 `CMCC Sync` 目录已不存在。当前配置：
  - 原始数据 `DATASET_PATH` = `E:\Medical Data\Neurosurgery\NF-PitNET`
  - 临时/派生 `DATASET_PATH_TEMP` = `F:\Data Temp\NF-PitNET`
  - 全部路径经 `scripts/_paths.py`（Python）/ `scripts/paths.R`（R）从项目 `.env` 读取，**脚本内无硬编码**。
- 软件版本与随机种子（Supplementary Methods 已写）：minfi v1.48.0 / R 4.3.3；scikit-learn 1.5.1，`random_state=42`；置换 20,000 次 seed=0；NN-mixing 置换 5,000 次。
- 数据来源、DOI、下载链接、本地布局：`docs/data-availability.md`。
- **原始数据一律不入库**（第三方版权），仓库只放代码 + 处理结果。
- 代码仓库：`github.com/PaperValidation/NF-PitNET-Methylation-Validation`（**已是 public**，2026-10-18 核实；稿件 Data Availability 中的链接可直接使用）。

---

## 10. 复核结论（2026-10-18）

- ✅ 数据完整性：原始数据 5 子目录 + 8 矩阵文件齐全；临时区 7 子目录完好；9 个关键数据文件经加载器全部可达。
- ✅ 产物完整性：11 个 summary JSON、17 个图件、全部补充表 CSV 在位。
- ✅ 代码可运行：`scripts/*.py` 全部 `py_compile` 通过。
- ✅ **数字级复核**：稿件全部关键断言逐项回溯源文件，**100% 吻合**（详见本文件第 2 节）。
- ✅ 稿件格式：6 项强制模块全部达标，16 条参考文献全部真实且被引用。
- ⚠️ 唯一未完成：6 处 `[AUTHOR_ACTION]`（需真人信息）。

**分析侧结论：可以进入成文/投稿包装阶段。**

---

## 11. 交接口径声明

本文件由 Biomni 生物医学专家团（分析侧）出具，**仅覆盖分析产物、证据链与诚实边界**。

按团队职责边界，以下工作**不在本团范围**，需由写作类专家团承接：
- 手稿正文撰写、结构重写、语言润色
- 摘要/关键词的重新撰拟（现有版本已合规，如需改写属写作范畴）
- Cover letter 定稿
- 审稿意见逐条回复的成文
- 期刊格式排版与投稿系统提交

若写作过程中**需要新的分析或数字**（如补充统计、新增亚组、重算某指标），请回到分析侧重新出数，**不要自行推算或改写既有数字**。
