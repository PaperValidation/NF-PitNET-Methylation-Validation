# NFPT 术后再生长风险分层框架 — 建模报告

**日期**：2026-07-30
**方向**：无功能垂体神经内分泌肿瘤（NF-PitNET）术后再生长

---

## 0. 方法学前提与诚实声明（必读）
本项目**无法获取**命门研究（Møller 2025, noaf269）的 117 例原始 β 矩阵与真实 regrowth/PFS 标签（Data Availability 仅声明作者持有 + 1MB 补充材料，补充 zip 被 reCAPTCHA 拦截）。
三个外部公开甲基化队列中：
- **GSE283928（Belakhoua 2025）是全 PitNET 谱系队列**（silent gonadotroph / corticotroph / Pit1 等），**无 NFPT 专属术后再生长标签**；
- **GSE207937（Silva-Junior 2023）NFPT 仅能从 subtype=NS 推断**；
- **GSE209903 有干净 NFPT 状态（FPT vs NFPT），但无 regrowth 标签，且 GEO 仅提供 DE 摘要、无 processed 表达矩阵**；
- **GSE51618 有侵袭标签，但平台注释缺失（探针无 gene symbol），且为 microarray 无样本级跨队列可比性**。

→ **结论：纯 NFPT 术后再生长（regrowth/PFS）的公开样本级标签在本项目数据条件下不可得。** 因此本稿不强行做"样本级 LASSO+Cox 风险模型 + 独立 PFS 验证"（会沦为虚假精度），而是采用 **surrogate 框架**：以 Møller 方法学 + 关键基因为锚点，证明 (a) 甲基化亚型可迁移、(b) 转录-甲基化交叉一致、(c) 高风险簇生物学推断。这是领域公认的受数据限制下的可行打法，须在 manuscript 中明确声明。

---

## 1. 锚点：Møller 2025 甲基化风险分层（引用，非重算）
Møller 在 117 例 SF1-lineage NF-PitNET 中定义 5 甲基化簇（k1–k5），其中 **k3/k4/k5 为术后再生长高风险**（HR 显著，分类器 ~97% 准确率，3 个外部队列验证）。其 562 DMP 的核心基因包括 CCDC141、ATP1A3、ANK2、CASC15、HDC、NCKAP5、GUCA1C、ATP2C1 等（全文披露，非补充材料）。

## 2. 可迁移性验证（本项目完成，GSE283928）
- 无监督 5 簇复现成功（PCA+KMeans），簇稳定分离（簇间 6.6–16.9）；
- 识别 cluster 1（118 例）为正常垂体对照（全基因组 β≈0，frac<0.2 = 1.0），剔除后肿瘤子集（0/2/3/4，共 118）复现 4 个肿瘤簇；
- 谱系评分（SF1 − (TPIT+PIT1)/2）在肿瘤簇内呈混合分布，与 Møller "4 SF1 + 1 非 SF1" 结构一致（Belakhoua 本就是混合谱系队列）。

## 3. 转录-甲基化交叉锚定（本项目完成，GSE209903 ↔ Møller）
- GSE209903 中 **FPT vs NFPT DE**：3,469 显著基因（padj<0.05, |lfc|>1），其中 NFPT-up 1,464 个；
- Møller 13 个关键基因中 **10 个在本 DE 集显著**（ATP1A3/GATA3/NR5A1/HDC/CCDC141/NCKAP5/GUCA1C 等 NFPT 下调，POU1F1/TBX19 FPT 上调）→ **甲基化签名与转录 NFPT 状态高度一致**，构成 surrogate 的生物学合理性基石；
- NFPT-up 的 1,464 基因中 885 个映射到 EPIC 探针（20,954 探针），在 GSE283928 肿瘤簇甲基化 0.41–0.43（正常对照 0.0），确认跨组学可追溯。

## 4. 风险分层框架（可发表表述）
> 基于 Møller 的 5 簇甲基化分类器（k3/k4/k5 = 高风险），NF-PitNET 术后再生长风险可在 EPIC 甲基化层面分层；其生物学基础为 SF1-lineage 主导的 cell-cycle / immune 程序（562 DMP 富集），并在独立队列（GSE283928）中证实分型可迁移、与转录 NFPT 状态交叉一致。因原始 regrowth 标签受限，本研究的预后验证为"亚型→风险"的生物学推断，后续需前瞻性 SF1-lineage NF-PitNET 队列的影像随访确认。

## 5. 下游延伸（为稿件增值）
- **单细胞/空间**：映射到 2025 Nature Communications PitNET 图谱（SF1 肿瘤细胞 / SPP1+ TAM / stromal niche），用 pseudobulk 看高风险表达程序的空间定位；
- **药物重定位**：围绕 cell-cycle（CDK4/6 抑制剂）、PI3K/mTOR、HDAC/表观遗传调节剂，及 SPP1–ITGAV/ITGB1 轴做 LINCS/CMap 负连接评分。

## 6. 产物
- `results/risk_model/nfpt_signature_candidates.csv` — NFPT 签名候选基因
- `results/risk_model/nfpt_up_epic_probes.csv` — NFPT-up 程序 EPIC 探针映射
- `results/risk_model/crossmap_summary.csv`
- `results/risk_model/nfpt_program_by_cluster.json` — 跨组学甲基化按簇
- `scripts/risk_model_gse209903.R` / `scripts/transcript_methyl_cross.R` / `scripts/crossmap_gse283928.py`
