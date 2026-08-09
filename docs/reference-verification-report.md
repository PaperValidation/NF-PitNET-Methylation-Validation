# 文献引用验证报告 (Reference Verification Report)

**稿件**: `docs/manuscript-v2.md` (Neuro-Oncology, Vancouver 编号制)
**验证日期**: 2026-08-07 · **收尾更新**: 2026-08-09
**最终文献数**: **16 条** (初稿 18 条; [2]、[4] 经复核为 AI 编造引用, 已于 2026-08-09 删除并重排编号)
**验证方法**: PubMed eutils (PMID + 题名检索) + Crossref API (题名/元数据), 经 SOCKS5 代理并行查询; 字段级比对作者首作者/题名/期刊/年/卷/期/页。
**说明**: 用户要求"作者信息暂时空着"——指稿件**署名作者栏 (byline)** 留空, 不影响本文对**参考文献列表**的核查。

---

## 结论速览 (2026-08-07 原始审计, 共 18 条)

- ✅ Verified (核心字段正确): **1** 条 — [7]
- ⚠️ Check suggested (核心正确但缺卷期页/题名缩写/年歧义): **13** 条 — [1],[3],[6],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18]
- ❌ Needs fix (内容错误, 必须修正): **4** 条 — [2],[4],[5],[8]

> 注: 初稿 18 条。其中 [2]、[4] 经 2026-08-09 二次复核为 **AI 编造引用 (不存在于数据库)**, 已从稿件删除; [5]、[8] 已按正确值修正。最终稿件 **16 条**, 全部真实且字段完整 (见"收尾状态"与下方修正后段落)。

---

## 收尾状态 (Resolution — 2026-08-09)

经二次复核, **原 [2] 与 [4] 并非"题名写错的真论文", 而是不存在于 PubMed/Crossref 的编造引用**。确认依据:

- **原 [2]** "Ntali G, Grossman A. … Pituitary. 2014" — Ntali 与 Grossman 共同署名的真实论文为 2013/2016 死亡率研究, **无 2014 年 Pituitary 的 NFPA 综述**; 且正文从未引用该编号 (孤儿)。→ **删除**。
- **原 [4]** "Sievers P et al. … Acta Neuropathol. 2020" — "Sievers + pituitary" 检索为空, Sievers 已索引论文均为 2024–2025 年垂体疾病的神经/认知方面, **无 methylation classifier**; 2020 年 PitNET 甲基化分型的唯一真实论文即 Neou 2020 (原 [9], 已存在于列表)。→ **删除**; 正文第 207 行对 "PitNET methylation classes" 的引用改指向 **Neou et al. 2020 (现 [7])** 与 Capper 2018 (现 [2])。

**处置动作** (commit `29bfef0`, 已推送 GitHub `main`):
1. 删除原 [2]、[4] 两条编造文献; 剩余 16 条重排为连续编号 **1–16** (双射映射, 语义无错位)。
2. 正文第 207 行: 引用 `(4)` → `(7)`; 作者 "Sievers et al. (2020)" → "**Neou et al. (2020)**"。
3. 校验: 16 条全部真实且被正文引用、编号连续、补充材料编号列表完好、无 `[VERIFY]` 残留。

> **经验教训**: ARS 写作技能自动生成的参考文献会编造 (hallucinate)。今后凡 generated 稿件, 必须逐条用 PubMed/Crossref 按"作者 + 期刊 + 年"核验**存在性与题名**, 不能仅留占位待用户补。

---

## ❌ 必须修正 (Critical)

| # | 问题 | 稿件当前值 | 正确值 (来源) | 严重度 |
|---|------|-----------|--------------|--------|
| **[5]** | 期号 + 页码全错 | *Neuro-Oncology* 2021;**23(3):441-453** | *Neuro Oncol* 2021;**23(8):1292-1303** (PMID 33631002, Crossref 一致) | 🔴 页码差 ≥5 + 期号错 |
| **[8]** | 作者/期刊/题名全错 (张冠李戴) | Silva-Junior **HM**, *Integrated methylome and transcriptome of PitNETs*, **Clin Epigenetics**, 2023 (GSE207937/GSE209903) | **Silva-Júnior RMP**, Bueno AC, Martins CS, et al. *Integrating Methylome and Transcriptome Signatures Expands the Molecular Classification of Pituitary Neuroendocrine Tumors.* **J Clin Endocrinol Metab.** 2022;108(6):1452-1463 (PMID 36504388; DOI 10.1210/clinem/dgac703) | 🔴 作者首作者错 + 期刊错 + 题名错 |
| **原 [2]** | 编造引用 (文献不存在) | Ntali G, Grossman A. *Clinical review: non-functioning pituitary adenomas.* *Pituitary.* **2014** | **已删除 (2026-08-09)** — 不存在于 PubMed/Crossref; Ntali/Grossman 真实论文为 2013/2016 死亡率研究; 正文亦未引用该编号 | 🔴 虚构引用, 删除 |
| **原 [4]** | 编造引用 (文献不存在) | Sievers P, et al. *DNA methylation-based classification of pituitary neuroendocrine tumors.* *Acta Neuropathol.* **2020** | **已删除 (2026-08-09)** — 不存在; "Sievers + pituitary" 检索为空, 其已索引论文均非 methylation classifier; 该论断真实出处为 Neou 2020 (现 [7]) 与 Capper 2018 (现 [2]) | 🔴 虚构引用, 删除 |

---

## ⚠️ 建议核对 (Warning — 核心正确, 需补元数据)

| # | 问题 | 当前值 → 建议修正 |
|---|------|------------------|
| [1] | 缺卷期页 + 题名过简 | *Neuro-Oncology* 2023 → **Neuro Oncol. 2023;25(12 Suppl 2):iv1-iv99** (PMID 37793125); 题名补全为 "CBTRUS Statistical Report: Primary Brain and Other Central Nervous System Tumors Diagnosed in the United States in 2016-2020" |
| [3] | 缺卷期页 (题名缩写可接受) | *Nature* 2018 → **Nature. 2018;555(7697):469-474** (PMID 29539639); 题名可补 "central nervous system tumours" |
| [6] | 首作者缩写错 + 缺卷期页 + 年歧义 | "Moller **AH**" → **Møller MW** (Morten Winkler Møller, PubMed/Crossref 一致); 补 **28(3):728-740**; 年: Crossref/DOI=2025 (稿件一致), PubMed epub=2026-03-01, 请确认最终卷年 |
| [9] | 缺卷期页 | 补 **Cancer Cell. 2020;37(1):123-134.e5** (PMID 31883967) |
| [10] | 缺卷期页 | 补 **Mol Cell Endocrinol. 2018;473:194-204** (PMID 29410024) |
| [11] | 缺卷期页 | 补 **Cancers (Basel). 2019;11(11):1650** (PMID 31731486) |
| [12] | 题名截断 + 缺卷期页 | "A next generation connectivity map" → 全名 "**A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles**"; 补 **Cell. 2017;171(6):1437-1452** (PMID 29195078) |
| [13] | 题名过简 + 缺卷期页 | "Enrichr" → 全名 "**Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool**"; 补 **BMC Bioinformatics. 2013;14:128** (PMID 23586463) |
| [14] | 缺卷期页 | 补 **Acta Neuropathol Commun. 2024;12(1):61** (PMID 38637883) |
| [15] | 缺卷期页 | 补 **Pituitary. 2023;26(4):474-481** (PMID 37428396) |
| [16] | 首作者缩写 + 缺卷期页 | "Shi **L**" → PubMed 为 **Shi C** (Chengzhang Shi); 补 **Neuro Oncol. 2020;22(8):1114-1125** (PMID 32246150) |
| [17] | 缺卷期页 | 补 **Oncotarget. 2016;7(43):69149-69158** (PMID 27708250) |
| [18] | 缺卷期页 | 补 **Neurosurg Focus. 2018;44(6):E9** (PMID 29852762) |

---

## ✅ 已验证无误

- **[7]** Belakhoua S, … Snuderl M. *DNA methylation profiling of pituitary neuroendocrine tumors identifies distinct clinical and pathological subtypes based on epigenetic differentiation.* Neuro Oncol. 2025;27(9):2341-2354. (PMID 40295206) — 首作者/题名/期刊/年/卷/期/页 全部一致, 完整正确。

> [9]–[18] (除上表标注的首作者缩写/缺卷期页外) 经 PMID 核验, 所引论文真实存在且与稿件指向一致, 仅投稿前需补齐卷期页。

---

## 修正后的参考文献段落 (与稿件 §6 一致 — 16 条终稿, 2026-08-09 收尾)

> 原 [2]/[4] 已删除, 编号重排为 1–16; 下表即稿件当前 §6 的最终内容。

```text
1. Ostrom QT, et al. CBTRUS Statistical Report: Primary Brain and Other Central Nervous System Tumors Diagnosed in the United States in 2016-2020. Neuro-Oncology. 2023;25(12 Suppl 2):iv1-iv99. (PMID 37793125)
2. Capper D, et al. DNA methylation-based classification of central nervous system tumours. Nature. 2018;555(7697):469-474. (PMID 29539639)
3. Mosella MS, Sabedot TS, Silva TC, et al. DNA methylation-based signatures classify sporadic pituitary tumors according to clinicopathological features. Neuro-Oncology. 2021;23(8):1292-1303. (PMID 33631002)
4. Møller MW, et al. DNA methylation profiling predicts postsurgical regrowth in SF1-lineage nonfunctioning pituitary neuroendocrine tumors. Neuro-Oncology. 2025;28(3):728-740. (PMID 41369112)
5. Belakhoua S, Vasudevaraja V, Schroff C, Galbraith K, Movahed-Ezazi M, Serrano J, Yang Y, Orringer D, Golfinos JG, Sen C, Pacione D, Agrawal N, Snuderl M. DNA methylation profiling of pituitary neuroendocrine tumors identifies distinct clinical and pathological subtypes based on epigenetic differentiation. Neuro-Oncol. 2025;27(9):2341-2354. (PMID 40295206)
6. Silva-Júnior RMP, Bueno AC, Martins CS, et al. Integrating Methylome and Transcriptome Signatures Expands the Molecular Classification of Pituitary Neuroendocrine Tumors. J Clin Endocrinol Metab. 2023;108(6):1452-1463. (PMID 36504388; GSE207937 / GSE209903)
7. Neou M, et al. Pangenomic Classification of Pituitary Neuroendocrine Tumors. Cancer Cell. 2020;37(1):123-134.e5. (PMID 31883967; E-MTAB-7762)
8. Kober P, et al. DNA methylation profiling in nonfunctioning pituitary adenomas. Mol Cell Endocrinol. 2018;473:194-204. (PMID 29410024; GSE115783)
9. Kober P, et al. The Role of Aberrant DNA Methylation in Misregulation of Gene Expression in Gonadotroph Nonfunctioning Pituitary Tumors. Cancers (Basel). 2019;11(11):1650. (PMID 31731486; GSE115783 follow-up)
10. Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell. 2017;171(6):1437-1452. (PMID 29195078; LINCS)
11. Chen EY, et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013;14:128. (PMID 23586463)
12. Cheng Y, et al. Therapeutic potential of targeting Nrf2 by panobinostat in pituitary neuroendocrine tumors. Acta Neuropathol Commun. 2024;12(1):61. (PMID 38637883; HDAC/panobinostat in PitNET)
13. Lin AL, et al. The treatment of aggressive prolactinomas with everolimus. Pituitary. 2023;26(4):474-481. (PMID 37428396; mTOR/everolimus in PitNET)
14. Shi C, et al. BRD4 as a therapeutic target for nonfunctioning and growth hormone pituitary adenoma. Neuro-Oncol. 2020;22(8):1114-1125. (PMID 32246150; BET/BRD4 in NFPA)
15. Zhang D, et al. Targeting the ERK pathway for the treatment of Cushing's disease. Oncotarget. 2016;7(43):69149-69158. (PMID 27708250; MEK, corticotroph-only)
16. Anderson E, et al. Regression of a nonfunctioning pituitary macroadenoma on the CDK4/6 inhibitor palbociclib: case report. Neurosurg Focus. 2018;44(6):E9. (PMID 29852762; CDK4/6, single NFPA case)
```

---

## 数据来源与方法备注

- 11 条带 PMID 的文献以 PubMed eutils `esummary` 核验 (题名/作者/期刊/年/卷/期/页)。
- 全部 18 条以 Crossref `query.bibliographic` 做二次比对; 题名相似度 (Jaccard) 用于排除张冠李戴。
- [2]/[4]/[8] 经题名精确检索、作者检索、期刊+年检索、宽泛主题检索共 5+ 轮交叉验证, 确认 [2]/[4] 查无、[8] 指向另一真实论文。
- 2026-08-09 二次复核确认 [2]/[4] 为编造引用, 已从稿件删除并重排编号 (见"收尾状态"); 正文误引 Sievers 2020 处改指 Neou 2020。
- 代理: 本机 SOCKS5 (127.0.0.1:8428); 临时验证脚本 `scripts/_ref_verify*.py` 已清理。
