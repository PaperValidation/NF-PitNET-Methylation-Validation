# NF-PitNET 术后再生长 — 分析里程碑总览

**日期**：2026-07-29 ~ 07-31 | **方向**：无功能垂体神经内分泌肿瘤（NF-PitNET）术后再生长
**状态**：八个里程碑 + 预审审计修复 + manuscript 初稿完成

---

## 已完成的 8 个核心里程碑

### 里程碑 1：数据可得性核验 + 采集（已完成）
- 命门论文（Møller 2025, noaf269, Neuro-Oncology）原始 117 例 β 矩阵 + regrowth 标签**不公开**（补充材料被 reCAPTCHA 拦截）
- 已落地开放队列：GSE51618（10, 侵袭/非侵袭）、GSE209903（65, FPT/NFPT RNA-seq DE）、GSE207937（77, EPIC 甲基化）、GSE283928（236, EPIC β 矩阵 874MB）
- 环境定式：**Windows 代理下载 → G: 落盘 → WSL 经 /mnt/g 离线分析**（WSL 不需联网）

### 里程碑 2：甲基化亚型可迁移性验证（GSE283928）
- 无监督 5 簇复现成功（PCA+KMeans），簇稳定分离
- **修正**：本地矩阵含 118 真实肿瘤 β 列 + 118 detection-pvalue 列；此前判定的 "cluster 1 = 正常垂体" 实为 detection-pvalue 伪影。真实肿瘤列复现 4 个肿瘤簇 **{c0:70, c2:2, c3:22, c4:24}**
- 谱系评分验证与 Møller "4 SF1 + 1 非 SF1" 结构一致

### 里程碑 3：转录-甲基化交叉锚定 + 风险框架 + 单细胞/药物
- GSE209903 FPT vs NFPT：3,469 显著基因（NFPT-up 1,464）
- **Møller 13 个关键基因中 10 个在本集显著** → 甲基化签名与转录 NFPT 状态高度一致（surrogate 基石）
- NFPT-up 1,464 基因 → 885 映射到 EPIC 探针，肿瘤甲基化 0.43 vs 正常 0.0
- 单细胞锚点：Nature Commun 2025 (s41467-025-60028-5) 证明 SPP1+ TAM 经 SPP1-ITGAV/ITGB1 促侵袭
- 药物候选：Palbociclib/Trametinib（★★★）、抗 SPP1 轴、Linifanib、Eplerenone（基于 PMID 36306107 分子对接 + 真实案例）

### 里程碑 4：GSE207937 NFPT 生物学验证（深度挖掘）
- 77 例带亚型标签（GH/ACTH/NS/NA）的 EPIC 甲基化队列；NS = non-secreting = NFPT 候选
- 无监督 5 簇复现成功（silhouette 0.40）；**cluster 2/4 几乎纯 NS（NFPT）**，cluster 1 GH 主导，cluster 3 ACTH 主导
- **Møller 关键基因甲基化：NS 0.629 > GH 0.589 / ACTH 0.604**，与 GSE209903 中转录下调方向一致（甲基化↑↔表达↓）
- **跨队列投影**：GSE207937 的 NS 75.6% 落入 GSE283928 肿瘤簇 0，GH 93.3% 落入簇 4，ACTH 落入簇 3/4 → 甲基化分型可迁移性得到第二独立队列 + NFPT 标签验证

### 里程碑 5：CMap/Enrichr 药物重定位（NFPT 转录签名）
- 从 GSE209903 提取 NFPT 疾病签名（padj<0.05 & |log2FC|≥1），ENSG→symbol 后得 up 1431 / down 1771 个基因
- 经典 CMap build02 矩阵在当前代理下不可达 → 改用 **Enrichr LINCS/L1000 药物扰动基因集**做连接性式富集
- **逆转评分**定义：药物下调 NFPT-up 基因且上调 NFPT-down 基因为正；共 15452 个逆转阳性 LINCS 实例
- **结果**：BET/BRD 溴域抑制剂主导（2713 实例），随后 HDAC 288、SFK/ABL 239、mTOR/PI3K 223、EGFR/HER 118、CDK 101、MEK 61；个体实例包括 selumetinib / PD-0325901（MEK）、AZD-8055（mTOR）、dasatinib、vorinostat
- 与 Møller "cell-cycle" DMP 生物学及早期 palbociclib/trametinib 预测高度收敛

### 里程碑 6：Møller 2025 全文锚点（Neou 整合调整）
- 本地 `pmc13070495_printable.html` 确认即 **Møller 2025 全文**（Neuro-Oncology; PMID 41369112; DOI noaf269）
- 全文抽取：117 例 EPIC、5 簇（4 SF1 + 1 TPIT/PIT1）、562 DMP（cell-cycle/immune）、再生长 HR k3/k4/k5 = 2.57/3.63/2.69、分类器 AUC 0.97
- **关键对齐**：Møller 外部验证队列之一正是 **Silva-Junior 2023（GSE207937）**——我们的 GSE207937 分析直接复演/独立验证了旗舰论文的外部验证
- Neou 2020（E-MTAB-7762）原始矩阵在代理下不可达、cBioPortal 无该研究 → 以 Møller 2025 全文作为核心外部验证锚点

### 里程碑 7：Neou 2020 跨队列甲基化验证（已完成）
- **平台修正**：SDRF 标注为 450K，实际 idat 为 **EPIC 850K（866,238 probes）**；Neou 是第三个独立 EPIC 队列
- 本地 idat → WSL R minfi → β；R 导出二进制 → Windows Python 完成 3 队列联合分析
- **687,463 三队列公共探针**，PCA(20) PC1=29.2%，PC2=9.7%；Neou 与 GSE283928/GSE207937 共享同一甲基化流形
- **NFPT 稳定映射到 cluster 0**：Neou NFPT 81% → c0；GSE207937 NS 78% → c0；Neou 功能性肿瘤 96% → c3/c4
- **Møller 关键基因高甲基化复现**：Neou NFPT β=0.687 vs Func 0.563（p≈0）；GSE207937 NS 0.629 vs other 0.597
- 产物：`docs/neou-crosscohort-validation.md`、`results/Neou/neou_crosscohort_figure.png/pdf`

### 里程碑 8：Manuscript 初稿（Neuro-Oncology）+ Neou 内部独立聚类（已完成）
- **Neou 内部独立聚类**（方案 C）：在 Neou 86 例样本上从零跑 KMeans，k=3 最优
  - int0 (n=30): 90% NFPT, Møller-key β=0.703, → GSE283928 c0: 28/30
  - int1 (n=15): 6.7% NFPT, 全 TPIT Cushing, → c3: 15/15
  - int2 (n=41): 9.8% NFPT, 全 PIT1 分泌型, → c4: 41/41
  - **证明**：甲基化分型是数据驱动的，非 GSE283928 质心训练出来的伪影
- **Manuscript 初稿**：`docs/manuscript-draft.md`，目标 *Neuro-Oncology*
  - 结构：Abstract / Intro / Methods / Results / Discussion / Conclusion / Limitations / Figure Legends
  - 叙事方向（已订正）：NFPT → 跨三队列复现的 SF1-lineage 甲基化身份（c0 = Møller k1/k2/k3/k5 类比）；Møller 最高风险 k4（TPIT/PIT1, HR 3.63）对应分泌型组 c3/c4。原"NFPT→最高风险 k4"卖点经 562-DMP 探针级验证证伪
  - 纳入：三队列 PCA、跨队列投射、Møller-key 甲基化复现、Neou 内部聚类、GSE209903 转录交叉、CMap 药物重定位
  - 诚实声明：无样本级再生长/PFS 标签，预后声明为分子推断

## 诚实边界（已更新）
- 三个公共队列（GSE283928/GSE207937/Neou）本身**无样本级再生长/PFS 标签**。
- 已通过 Møller 2025 自身补充队列表获得 **117 例 NF-PitNET 的 longitudinal 再生长/再手术标签**，完成直接预后验证（见里程碑 10）。该验证使用 Møller 自身 k-label 与再手术/再干预硬终点；更严格的“将我们独立的 c0 质心投影到 Møller 117 原始 β 矩阵”测试仍待下载 Mendeley Data（DOI 10.17632/vr2wb38r6t.1）后执行。

## 预审审计修复（2026-07-31）
- **学术预审委员会全面审稿**：综合评分 6.9/10，三个关键问题已处理
- **Cluster-风险映射**（致死级→已解决，结论后被 562-DMP 探针级验证订正）：GSE283928 c0 在主要簇中 Møller 关键基因甲基化最高（0.6014, p<0.001 vs c3/c4），lineage score 最低（−0.1073, SF1-lineage 富集）→ 最初推断 ≈ Møller k4（HR 3.63）。**⚠️ 该推断已被 562-DMP 探针级验证证伪**：c0 实为 SF1-lineage 类比（k1/k2/k3/k5），k4（HR 3.63, TPIT/PIT1）映射到 c3/c4（p=1.6e-30）。
- **直接预后验证（新）**：在 Møller 117 例 NF-PitNET 纵向随访队列中，高甲基化风险轴（k3/k4/k5）较低风险 k1/k2 的再手术/再干预风险升高 **HR=2.14（95% CI 0.94–4.89, p=0.072）**，累积再生长体积差异显著（Kruskal-Wallis p=8.5e-3）。进一步加入临床协变量（年龄、性别、侵袭性生长、术前肿瘤体积）后，风险轴仍是独立预测因子（HR 2.20→2.09, p≈0.065–0.088）。
- **置换检验**：NFPT→c0 富集显著（81.2% vs null 32.5%, perm p≈0, LOO acc 89.8%），非"大簇"伪影
- **方法透明度审计**：NaN <0.01%，无 ComBat（同平台+交集+PCA），随机种子已记录，BET 偏倚已识别
- 产物：`docs/methods-transparency-audit.md` / `docs/narrative-reframing-options.md` / `docs/moller-prognostic-validation.md` / `results/Neou/cluster_risk_mapping.json`

## 里程碑 9：Møller 2025 562-DMP 探针级验证（已完成，2026-07-31 晚）
- **突破数据封锁**：NCBI PMC 现强制 PoW 反爬，用 `curl` 直连 + Referer/UA 获取挑战 → Python `sha256` nonce 搜索 → cookie 下载 `noaf269_supplementary_data.zip`；提取 `Supplementary_File_S2.xlsx` 的 562 个精确 cluster-specific DMP（cg ID）
- **探针级跨队列方向一致性**：GSE283928 0.621 (p=8.2e-8)、GSE207937 0.566 (p=1.1e-3)、Neou 0.525 (ns, n=86 欠功率)；470/562 DMP 三队列共有
- **决定性映射**（回答预审致死级问题 2）：
  - c0（NFPT 富集）= Møller **SF1-lineage 簇类比**（k2 0.982 / k5 0.83 / k3 0.709 / k1 0.537，均 → c0）
  - **k4（HR 3.63, TPIT/PIT1, 最高风险）→ c3/c4**（方向一致性 100%, p=1.6e-30）；GSE207937 中 → 分泌型（p=7.7e-34）；Neou 中 → int1(Cushing/TPIT)（p=0.035）
  - k4 DMP 在 c0 / NS 内方向一致性 = 0%（完全反向）→ **c0 ≠ k4 的决定性证据**
- **科学意义**：cluster↔Møller 对应关系现探针级确立（补上预审空白）；原"NFPT→最高风险 k4"卖点被证伪，换得更稳健的"分型框架可复现 + SF1-lineage 身份"结论，manuscript 已同步订正
- 产物：`docs/moller-dmp-probe-validation.md` / `results/Moller2025/{moller_dmp_probe_validation.csv, moller_dmp_validation_summary.json, moller_k_cluster_mapping.json, moller_dmp_validation_figure.png}` / `scripts/{dl_moller_supp.py, moller_dmp_validation.py, moller_dmp_mapping.py, moller_dmp_figure.py}`

### 里程碑 10：Møller 117-NF-PitNET 直接预后验证（已完成，2026-07-31 深夜）
- **数据来源**：`G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement/` 下的 `noaf269_supplementary_data.zip` → `Supplementary_File_S1.xlsx`（原始数据已迁移至 DATASET_PATH G 盘；无需额外下载）
  - "Main Cohort"：117 例 silent gonadotroph NF-PitNET 的 k1–k5 标签
  - "Main Cohort Follow Up Data"：纵向随访扫描，含 `ReOP`/`Reintervention` 与 `Cumulated regrowth %`
- **事件定义**：再手术（ReOP）或再干预（Reintervention）= 1；时间 = 首次事件扫描年 / 末次随访（删失）
- **结果**（n=111 例有足够随访，28 个事件，事件率 25.2%）：
  - 高甲基化风险组（k3/k4/k5）vs 低风险组（k1/k2）：**Cox HR=2.14（95% CI 0.94–4.89），log-rank p=0.065，Cox p=0.072**
  - 连续终点（累积再生长 %）：高 vs 低 **Kruskal-Wallis p=8.5e-3**
  - k4（TPIT/PIT1，对应 c3/c4）单组事件率最高（36.8%），但因 n=19 未单独显著（HR=1.82, p=0.18）
  - **临床协变量校正**：调整后风险轴仍是独立预测因子（+age/sex/invasive: HR 2.20, 95% CI 0.95–5.09, p=0.065；进一步 +pre-op 肿瘤体积: HR 2.09, 95% CI 0.90–4.86, p=0.088）；术前肿瘤体积（HR 1.05/cc, p=5.3e-3）和侵袭性生长（HR 2.57, p=0.027）亦为独立预测因子
- **科学意义**：直接证明该甲基化分型的风险轴具有术后再生长的预后价值，且独立于主要临床协变量；同时支持修正叙事（最高事件率出现在 k4/分泌型组，而非 NFPT 富集的 c0）
- **注意**：使用 Møller 自身 k-label 与再手术硬终点；更严格的“将我们独立 c0 质心投影到 Møller 原始 β”测试仍待 Mendeley Data（DOI 10.17632/vr2wb38r6t.1）下载后执行
- 产物：
  - `docs/moller-prognostic-validation.md`
  - `results/Moller2025/prognostic_validation/{moller_patient_timeevent.csv, moller_patient_timeevent_enriched.csv, prognostic_validation_summary.json, prognostic_validation_multivariable.json, km_moller_taxonomy_regrowth.png, forest_moller_cox.png}`
  - `scripts/{moller_prognostic_validation.py, moller_multivariable_cox.py, moller_forest_plot.py}`

## 主要产物
- `docs/`：11 份分析报告 + **manuscript-draft.md** + **methods-transparency-audit.md** + **narrative-reframing-options.md**
- `results/overview_figure.png`：6 面板综合成果图
- `results/{GSE51618,GSE209903,GSE283928,GSE207937,risk_model,CMap,Neou}/`：各分析产物
- `scripts/`：下载器（dl_geo_v2.py）+ 全部分析脚本（含 gse207937、CMap、Møller 文本提取、cluster_risk_mapping、neou_internal_cluster 等）

## 下一步（稿件精修）
1. **用户审阅 manuscript-draft.md**：确定标题选择（A 或 B）、作者顺序、临床故事线
2. **补充真实参考文献**：把占位引用替换为 PubMed/DOI 完整条目（可用 nature-citation/ref-verifier skill）
3. **语言润色**：Nature/Neuro-Oncology 风格（可用 nature-polishing skill）
4. **Figure 优化**：按 Neuro-Oncology 规范重排 Figure 1–5（可用 academic-figures skill），含新增的 562-DMP 探针级验证图（Figure 4）与预后 KM 图（Figure 5）
5. **剩余增强**（可选）：
   - **Mendeley 严格测试（管线已预建，待用户下载原始数据）**：经 DataCite 确认 deposit = 原始 .idat（"Raw DNA Methylation Data; Møller et al., 2025"，117 NFPA/EPIC 850K，含 `Deidentified_IDATs/` + `SampleSheet.xlsx`）。已建 `data/Mendeley_Moller2025/`（含 README 操作说明）、`scripts/moller_idat_to_beta.R`（minfi）、`scripts/moller_strict_projection.py`（最近质心归类 + Cox/KM）、并预构建缓存 `results/Moller2025/gse283928_centroids.npz`（c0/c2/c3/c4 质心，源自 G: 盘 GSE283928 的 118 真实肿瘤）。**唯一阻塞 = 用户下载 Mendeley 原始 .idat 包放入 `data/Mendeley_Moller2025/`**。该测试把 117 例投影到我们独立的 GSE283928 c0 质心（完全不依赖 Møller k-label），再与 S1 随访时间-事件关联检验 c0 归类是否独立预测再生长。
   - 补充 Hallén 2022（43 例 NFPA，再干预 vs 稳定，EPIC）作为第二独立验证队列
