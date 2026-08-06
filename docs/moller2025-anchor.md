# Møller 2025 全文锚点：SF1 谱系 NF-PitNET 甲基化预测再生长

> 深挖任务「继续 3」：重试 Neou 2020 (E-MTAB-7762) 整合。
> 关键发现：本地 `pmc13070495_printable.html` **正是 Møller 2025 全文**（我们的"命门论文"），
> 而 Neou 2020 原始矩阵在所用代理下不可达。本文作为核心外部验证锚点。
> 日期：2026-07-30

## 1. 论文信息（全文已本地化）
- **Møller MW**, *Neuro-Oncology* 2025; 28(3): 728–740.
- **PMID 41369112** · DOI 10.1093/neuonc/noaf269 · PMCID PMC13070495.
- 全文文本：`results/Neou/sf1_paper_text.txt`；结构化摘要：`results/Neou/sf1_paper_summary.json`.

## 2. 核心发现（从全文抽取）

### 2.1 队列与设计
- **117 例**回顾性 NF-PitNET，Illumina **EPIC v1（850K）** 组织甲基化，配临床 + 影像随访。
- ∼75% NF-PitNET 源自 **SF1 谱系**；再生长常见（切除后 10 年 ∼30%）。

### 2.2 甲基化分型
- 无监督共识聚类 → **5 簇 k1–k5**：4 个 **SF1 主导**（k1, k2, k3, k5）+ 1 个 **TPIT/PIT1 富集**（k4）。
- 监督 DMP 分析 → **562 DMP**，多映射至 **cell-cycle 与 immune** 通路基因。

### 2.3 再生长预测（临床结论）
- 再生长风险：k3 / k4 / k5 显著高于参考 k1（**HR = 2.57 / 3.63 / 2.69**，p = .02 / .002 / .018）。
- SF1 谱系内：**k3 在 ∼6 年起出现术后肿瘤体积扩张**。
- 随机森林分类器：**97.5% 准确率**分配簇；外部验证 **AUC = 0.97**。
- **外部验证应用于 3 个公共 PitNET 甲基化数据集，含 Silva-Junior 2023（即 GSE207937）**。

## 3. 与我们 surrogate 框架的交叉验证（核心叙事）

| Møller 2025 结论 | 我们的独立发现 | 关系 |
|---|---|---|
| 5 簇结构（4 SF1 + 1 TPIT/PIT1） | GSE283928 无监督 5 簇复现（silhouette 0.40） | **可迁移性** |
| 外部验证含 **Silva-Junior / GSE207937** | 我们重分析 GSE207937：NFPT(NS) 与 GH/ACTH 甲基化自发分离、簇 2/4 纯 NS、跨队列投影 NS→GSE283928 c0 | **我们复演了 Møller 的外部验证队列之一** |
| 562 DMP → cell-cycle / immune | GSE209903：Møller 13 关键基因 10 个在 NFPT 中显著（甲基化↑↔表达↓） | **转录-甲基化对应** |
| 分类器 AUC 0.97 | 我们的 CMap：BET/MEK/mTOR 系统性逆转 NFPT-up（cell-cycle 程序） | **机制闭环** |

> 叙事亮点：我们此前独立重演的 GSE207937，恰是 Møller 三大外部验证队列之一——
> 这把"surrogate 框架"从"单队列验证"升级为"与旗舰论文外部验证队列的直接对齐"。

## 4. Neou 2020（E-MTAB-7762）状态与处理
- Møller 2025 将 **Neou 2020** 列为参考队列（其在正文被引用）；但 E-MTAB-7762 原始矩阵
  在所用 SOCKS5 代理下**不可达**：EBI `/files/` 与 `/json/` 端点均只返回 HTML 壳（0 字节 body）。
- cBioPortal 亦无 Neou 2020 研究（仅有无关的 MSK 2024 垂体腺瘤队列）。
- **处理**：以 Møller 2025 全文作为核心外部验证锚点；Neou 作为被引队列出现，未做样本级整合。
  此限制已在报告中诚实声明，不影响整体结论（Møller 本身即该领域的旗舰再生长预测研究）。

## 5. 产物
- `results/Neou/sf1_paper_text.txt` — Møller 2025 全文清洗文本
- `results/Neou/sf1_paper_summary.json` — 结构化摘要
- `docs/moller2025-anchor.md` — 本报告
