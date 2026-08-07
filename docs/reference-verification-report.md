# 文献引用验证报告 (Reference Verification Report)

**稿件**: `docs/manuscript-v2.md` (Neuro-Oncology, Vancouver 编号制, 共 18 条文献)
**验证日期**: 2026-08-07
**验证方法**: PubMed eutils (PMID + 题名检索) + Crossref API (题名/元数据), 经 SOCKS5 代理并行查询; 字段级比对作者首作者/题名/期刊/年/卷/期/页。
**说明**: 用户要求"作者信息暂时空着"——指稿件**署名作者栏 (byline)** 留空, 不影响本文对**参考文献列表**的核查。

---

## 结论速览

- ✅ Verified (核心字段正确): **1** 条 — [7]
- ⚠️ Check suggested (核心正确但缺卷期页/题名缩写/年歧义): **13** 条 — [1],[3],[6],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18]
- ❌ Needs fix (内容错误, 必须修正): **4** 条 — [2],[4],[5],[8]

> 注: [9]–[18] 除缺卷期页外, 首作者/题名/期刊/年经 PMID 核验均正确; 卷期页在投稿终稿由期刊要求补齐, 列为 ⚠️ 提醒。

---

## ❌ 必须修正 (Critical)

| # | 问题 | 稿件当前值 | 正确值 (来源) | 严重度 |
|---|------|-----------|--------------|--------|
| **[5]** | 期号 + 页码全错 | *Neuro-Oncology* 2021;**23(3):441-453** | *Neuro Oncol* 2021;**23(8):1292-1303** (PMID 33631002, Crossref 一致) | 🔴 页码差 ≥5 + 期号错 |
| **[8]** | 作者/期刊/题名全错 (张冠李戴) | Silva-Junior **HM**, *Integrated methylome and transcriptome of PitNETs*, **Clin Epigenetics**, 2023 (GSE207937/GSE209903) | **Silva-Júnior RMP**, Bueno AC, Martins CS, et al. *Integrating Methylome and Transcriptome Signatures Expands the Molecular Classification of Pituitary Neuroendocrine Tumors.* **J Clin Endocrinol Metab.** 2022;108(6):1452-1463 (PMID 36504388; DOI 10.1210/clinem/dgac703) | 🔴 作者首作者错 + 期刊错 + 题名错 |
| **[2]** | 文献不存在 (PubMed/Crossref 均查无) | Ntali G, Grossman A. *Clinical review: non-functioning pituitary adenomas.* *Pituitary.* **2014** | **未找到** — 以 Ntali+Grossman+Pituitary+2014 及题名多种检索均无对应论文; 最接近的实为 Zatelli 2018 / Penn 2017 / Greenman 2015 等综述, 均非该署名 | 🔴 虚构/错引, 须替换 |
| **[4]** | 文献不存在 (Sievers 无此论文) | Sievers P, et al. *DNA methylation-based classification of pituitary neuroendocrine tumors.* *Acta Neuropathol.* **2020** | **未找到** — Sievers 已索引论文均为垂体疾病的精神/认知方面, 无 methylation classifier; 跨 PubMed/Crossref 多策略检索均无此标题论文 (2020 前后) | 🔴 虚构/错引, 须替换 |

> **[2] 与 [4] 的处理建议**: 这两条在数据库中均不存在, 不能自动纠正。请核实您真正想引用的论文:
> - 若 [4] 意在指"通用 CNS 甲基化分类器", 那已是 **[3] Capper 2018**, 建议删除 [4] 或改为正确的 PitNET 专用分类文献。
> - [2] 的非功能垂体腺瘤临床综述, 请确认正确出处 (例如 Greenman Y, Stern N. *Optimal management of non-functioning pituitary adenomas.* Endocrine. 2015;50(1):51-55, 或您实际参考的综述)。

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

## 修正后的参考文献段落 (可直接替换稿件 §6)

```text
1. Ostrom QT, Price M, Neff C, et al. CBTRUS Statistical Report: Primary Brain and Other Central Nervous System Tumors Diagnosed in the United States in 2016-2020. Neuro Oncol. 2023;25(12 Suppl 2):iv1-iv99.
2. [VERIFY — 原 "Ntali G, Grossman A. Clinical review: non-functioning pituitary adenomas. Pituitary. 2014" 在 PubMed/Crossref 查无; 请替换为正确综述]
3. Capper D, Jones DTW, Sill M, et al. DNA methylation-based classification of central nervous system tumours. Nature. 2018;555(7697):469-474.
4. [VERIFY — 原 "Sievers P, et al. DNA methylation-based classification of pituitary neuroendocrine tumors. Acta Neuropathol. 2020" 查无; 若指通用分类器见 [3], 否则替换]
5. Mosella MS, Sabedot TS, Silva TC, et al. DNA methylation-based signatures classify sporadic pituitary tumors according to clinicopathological features. Neuro Oncol. 2021;23(8):1292-1303.
6. Møller MW, Herrgott GA, Andersen MS, et al. DNA methylation profiling predicts postsurgical regrowth in SF1-lineage nonfunctioning pituitary neuroendocrine tumors. Neuro Oncol. 2025;28(3):728-740.
7. Belakhoua S, Vasudevaraja V, Schroff C, et al. DNA methylation profiling of pituitary neuroendocrine tumors identifies distinct clinical and pathological subtypes based on epigenetic differentiation. Neuro Oncol. 2025;27(9):2341-2354.
8. Silva-Júnior RMP, Bueno AC, Martins CS, et al. Integrating Methylome and Transcriptome Signatures Expands the Molecular Classification of Pituitary Neuroendocrine Tumors. J Clin Endocrinol Metab. 2022;108(6):1452-1463.
9. Neou M, Villa C, Armignacco R, et al. Pangenomic Classification of Pituitary Neuroendocrine Tumors. Cancer Cell. 2020;37(1):123-134.e5.
10. Kober P, Boresowicz J, Rusetska N, et al. DNA methylation profiling in nonfunctioning pituitary adenomas. Mol Cell Endocrinol. 2018;473:194-204.
11. Kober P, Boresowicz J, Rusetska N, et al. The Role of Aberrant DNA Methylation in Misregulation of Gene Expression in Gonadotroph Nonfunctioning Pituitary Tumors. Cancers (Basel). 2019;11(11):1650.
12. Subramanian A, Narayan R, Corsello SM, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell. 2017;171(6):1437-1452.
13. Chen EY, Tan CM, Kou Y, et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013;14:128.
14. Cheng Y, Dai Y, Tang H, et al. Therapeutic potential of targeting Nrf2 by panobinostat in pituitary neuroendocrine tumors. Acta Neuropathol Commun. 2024;12(1):61.
15. Lin AL, Geer EB, Lala N, et al. The treatment of aggressive prolactinomas with everolimus. Pituitary. 2023;26(4):474-481.
16. Shi C, Ye Z, Han J, et al. BRD4 as a therapeutic target for nonfunctioning and growth hormone pituitary adenoma. Neuro Oncol. 2020;22(8):1114-1125.
17. Zhang D, Bergsneider M, Wang MB, et al. Targeting the ERK pathway for the treatment of Cushing's disease. Oncotarget. 2016;7(43):69149-69158.
18. Anderson E, Heller RS, Lechan RM, et al. Regression of a nonfunctioning pituitary macroadenoma on the CDK4/6 inhibitor palbociclib: case report. Neurosurg Focus. 2018;44(6):E9.
```

---

## 数据来源与方法备注

- 11 条带 PMID 的文献以 PubMed eutils `esummary` 核验 (题名/作者/期刊/年/卷/期/页)。
- 全部 18 条以 Crossref `query.bibliographic` 做二次比对; 题名相似度 (Jaccard) 用于排除张冠李戴。
- [2]/[4]/[8] 经题名精确检索、作者检索、期刊+年检索、宽泛主题检索共 5+ 轮交叉验证, 确认 [2]/[4] 查无、[8] 指向另一真实论文。
- 代理: 本机 SOCKS5 (127.0.0.1:8428); 脚本见 `scripts/_ref_verify*.py` (临时验证脚本, 可删)。
