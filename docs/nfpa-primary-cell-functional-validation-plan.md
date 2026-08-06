# NFPA 原代细胞功能验证路线图（Translational Validation Plan）

> 关联：manuscript §3.9 / §4（CMap 药物重定位，已按 deep-research 结论诚实降级）
> 日期：2026-08-01
> 目标：把"SF1-lineage 甲基化身份 + CMap 可成药节点"的**计算假说**推进到**NFPA 原代细胞功能性证据**。

---

## 1. 我们要验证的假说

1. **谱系身份药敏假说（核心，目前零直接证据）**：NF-PitNET 的 SF1-lineage（c0）甲基化身份是否预测对 HDAC / mTOR / BET 类抑制剂的敏感性？
   - 现状：谱系身份只有**预后/分型**证据（Møller 2025 k3/k4/k5 再生长分层；本项目的 c0 映射）。**没有任何研究检验 SF1-lineage 状态是否预测对任何抑制剂的响应**——这是纯假说，manuscript 已明确不写预测性声明。
2. **候选药优先度假说**：panobinostat（HDAC）/ everolimus（mTOR）/ BET（JQ1）在 NFPA 原代细胞中应有可测的抗增殖/促凋亡效应；MEK / CDK4/6 仅作低置信探索。

---

## 2. 可用的原代 / 患者来源模型（文献已核实）

| 模型 | 来源 | 规模 / 特点 | 适用 |
|---|---|---|---|
| **患者来源 PitNET 原代细胞 panel** | Li *et al.* 2024 (PMID 38637883) | **9 例** patient-derived PitNET 原代细胞 + 高通量药筛，**panobinostat（HDAC）是最强单药类** | 直接可用面板，优先 |
| **NFPA 原代培养 + 异种移植** | Shi *et al.* 2020 (BET 在 NFPA) | NFPA 原代细胞 + 小鼠异种移植，测 BET（JQ1/OTX015）增殖/凋亡 | BET 假说验证 |
| **标准 NFPA 原代短期培养** | 通用垂体瘤原代协议 | 手术标本酶解 → 短期原代（数天~数周），做增殖/凋亡/激素 | 基础平台 |
| **GH3 / AtT-20 / MMQ 细胞系** | 文献常备 | 已建系，但非 NFPA 专属（GH3=GH, AtT-20=ACTH） | 机制对照，非 NFPA 证据 |
| **患者来源类器官 (PDO) / PDX** | **极稀缺，未见 NFPA 专属发表** | 模型瓶颈 | 中长期目标 |

> 模型瓶颈：PitNET（尤其 NFPA）患者来源类器官 / PDX 几乎未见发表；原代细胞寿命短、扩增难。这是本路线图的最大实操约束。

---

## 3. 候选药优先级（承接 manuscript 降级结论）

| Tier | 药物类（代表） | PitNET 已有证据 | NFPA 专属证据 | 本验证定位 |
|---|---|---|---|---|
| **1（高）** | HDAC（panobinostat） | Li 2024 原代面板最强单药；Zhang 2021 CUDC-907 (PI3K/HDAC) Cushing 异种移植 ↓65% | Li 2024 含 NFPA 来源 | **主攻** |
| **1（高）** | mTOR（everolimus） | Di Pasquale 2018 GH3/人 GH 培养；Lin 2023 侵袭性泌乳素瘤临床 | Zatelli 2024 综述提及 NFPA 培养探索 | 主攻 |
| **1（高）** | BET（JQ1/OTX015） | Shi 2020 NFPA 原代+xenograft；Lines 2020 AtT20 增殖↓95% | Shi 2020 直接 NFPA | 主攻 |
| **2（低置信）** | MEK（trametinib） | 仅 corticotroph（Zhang 2016 binimetinib） | **无** | 探索，不承诺 |
| **2（低置信）** | CDK4/6（palbociclib） | Anderson 2018 单例 NFPA 偶然退缩 | 仅 1 例个案 | 探索，不承诺 |

---

## 4. 实验设计（proposed）

**4.1 模型准备**
- 优先接入 Li 2024 式 9 例 patient-derived 原代 panel（或自建：手术 NFPA 标本 → 胶原酶消化 → 种植于基质胶/纤连蛋白，血清-free 或低血清短期培养）。
- 对每个供体**同步做 EPIC 甲基化分型**（用本项目 c0 质心 / 562-DMP 分类器），标记其 SF1-lineage（c0）身份 —— 这是检验"谱系身份→药敏"假说的关键衔接。

**4.2 处理与终点**
- 药物梯度：1 nM – 10 μM，DMSO 对照；72–120 h。
- 主要终点：
  - 增殖：WST-1 / EdU 掺入 / 活细胞计数
  - 凋亡：caspase-3 活性 / Annexin V
  - 细胞周期：PI 染色（HDAC/CDK4/6 应引起 G1 阻滞）
  - 激素分泌（若供体保留）：LH/FSH/E2（SF1-lineage 相关）
- 阳性对照：各药类已知敏感模型（如 BET 用 AtT-20）。

**4.3 假说检验矩阵**
- 主检验：Tier-1 药（panobinostat / everolimus / JQ1）在 NFPA 原代中的剂量-响应曲线（IC50 / 最大抑制率）。
- 分层检验（核心）：按供体 c0（SF1-lineage）身份分层，比较 HDAC/BET 敏感性差异 → 直接回答"谱系身份→药敏"假说（目前零证据，本实验首次检验）。

**4.4 统计**
- 每药每供体 ≥3 复孔；IC50 用 4-parameter logit；分层比较用配对/混合效应模型；多重比较校正。

---

## 5. 与我们分型的衔接（闭环价值）

- 用本项目已建立的 **c0 质心投影** 与 **562-DMP 分类器** 对每个原代细胞供体做甲基化分型，使"计算分型"与"功能药敏"在同一供体上配对 → 把 manuscript 的甲基化分型从 surrogate 直接升级为有功能意义的生物学实体。
- 若 c0/SF1-lineage 身份确实预测 HDAC/BET 敏感性，则 manuscript 的"SF1-lineage → 药敏"可从**假说**升级为**证据**，是重大卖点；若无关，则诚实限定为"谱系分型独立于药敏，仅预后相关"。

---

## 6. 缺口与风险

1. **SF1-lineage→药敏 零直接证据**：本路线图核心假说未经任何研究检验，可能阴性。
2. **模型瓶颈**：NFPA 原代细胞寿命短；PDO/PDX 稀缺，难以规模化药筛。
3. **NFPA 生长慢**：cytostatic 药（HDAC/CDK4/6）效应需长期 assay，短期读数可能低估。
4. **血脑屏障 / 垂体位置**：CNS 渗透药（尤其大分子）在垂体局部递送挑战，影响体内外一致性。
5. **跨平台缺口**：当前外部验证三队列均为 EPIC 850K（同平台）；HM450 队列 GSE115783 被 NCBI 阻断下载，跨平台验证待补。

---

## 7. 参考文献（deep-research 已核实）

- Li *et al.* 2024 — 9 例 patient-derived PitNET 原代细胞 HTS，panobinostat 最强单药类. PMID 38637883.
- Shi *et al.* 2020 — BET 抑制剂在 NFPA 原代 + 异种移植. (BET/NFPA 主攻证据)
- Zatelli *et al.* 2024 — NFPA 药物治疗证据综述（证据层级）. PMID 39736867.
- Di Pasquale *et al.* 2018 — everolimus 在 GH3 / 人 GH 分泌培养. (mTOR)
- Lin *et al.* 2023 — everolimus 侵袭性泌乳素瘤临床. (mTOR 临床)
- Anderson *et al.* 2018 — palbociclib 单例 NFPA 偶然退缩. (CDK4/6 低置信)
- Zhang *et al.* 2016 — binimetinib (MEK) corticotroph 模型. (MEK 仅 corticotroph)
- Zhang *et al.* 2021 — CUDC-907 (PI3K/HDAC) Cushing 异种移植. (HDAC 强证据)
- Sood *et al.* 2024 — SF1 IHC 对 gonadotroph 肿瘤敏感性/特异性. (谱系标志物)
- Møller *et al.* 2025 — 117 NF-PitNET 甲基化 5 簇 + 再生长 HR. PMID 41369112.
