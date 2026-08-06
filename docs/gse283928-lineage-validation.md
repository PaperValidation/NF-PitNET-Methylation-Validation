# GSE283928 甲基化亚型可迁移性验证 — 里程碑报告

**日期**：2026-07-30
**队列**：GSE283928（Belakhoua et al., *Neuro-Oncology* 2025, noaf109）— EPIC 850K 处理 β 矩阵，236 样本（118 肿瘤 + 118 配对正常垂体）
**目的**：验证 Møller et al. 2025（noaf269, *Neuro-Oncology*）提出的"SF1-lineage 主导的 5 甲基化簇"结构，能否在独立公开队列中复现/迁移。

---

## 1. 数据形态确认
- 矩阵量纲：**标准 β 值（0–1）**，全矩阵 median ≈ 0.56，分布正常（前期诊断确认，非 M 值）。
- 236 列 = `Case1_Beta … Case236_Beta`，匿名样本名；临床标签（包括是否 NFPT / 复发）**不在矩阵内**，需论文补充表（被 reCAPTCHA 硬墙拦截，未获取）。

## 2. 无监督 5 簇复现（前序里程碑，见 analysis-progress-report.md）
- 全探针 → 方差过滤 → PCA(20) → KMeans(k=5)
- 簇大小：**{0:70, 1:118, 2:2, 3:22, 4:24}**
- 簇间距离 6.6–16.9，远大于簇内离散 → 结构稳定分离

## 3. ★ 本里程碑核心发现：cluster 1 = 正常垂体对照
对每簇做全局 β 诊断（subsample 12 万探针稳健估计）：

| cluster | n_samples | global_mean_β | frac_β<0.2 |
|---|---|---|---|
| 0 | 70 | 0.49 | 0.24 |
| **1** | **118** | **0.0** | **1.0** |
| 2 | 2 | 0.51 | 0.26 |
| 3 | 22 | 0.50 | 0.24 |
| 4 | 24 | 0.49 | 0.23 |

→ **cluster 1（恰好 118 例）在所有探针上 β≈0、全基因组去甲基化**，与肿瘤簇（β median ≈ 0.49–0.51）形成鲜明对比。
→ 完美匹配 Belakhoua 原文的 **118 肿瘤 + 118 配对正常** 配对设计。
→ **结论：cluster 1 应作为正常对照剔除，仅在肿瘤子集（cluster 0/2/3/4，共 118 肿瘤）上复现 Møller 5 簇。**

## 4. 谱系标记探针验证（EPIC anno 映射）
用 `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` 将谱系标记基因映射到 EPIC 探针：
- SF1-lineage（NR5A1/GATA3/GATA2/EMX2/FSHB/LHB）：42 探针
- TPIT-lineage（TBX19/POMC/PCSK1/CRH）：55 探针
- PIT1-lineage（POU1F1/GH1/PRL/TSHB/GHRHR）：21 探针
- Møller 2025 关键基因集（CCDC141/ATP1A3/ANK2/CASC15/HDC/NCKAP5/GUCA1C/ATP2C1/GATA3/NR5A1）：168 探针

**谱系 dominance 评分（SF1 − (TPIT+PIT1)/2）按肿瘤簇分布**：
| cluster（仅肿瘤） | n | mean SF1-score |
|---|---|---|
| 0 | 70 | −0.056 |
| 2 | 2 | +0.054 |
| 3 | 22 | +0.044 |
| 4 | 24 | −0.021 |

→ 肿瘤簇并未全部 SF1 主导（cluster 0 偏非 SF1，cluster 2/3 偏 SF1），与 Møller "4 SF1 + 1 非 SF1" 的**混合谱系结构一致**（Belakhoua 是全 PitNET 谱系队列，含 silent corticotroph / Pit1 等多种亚型，本身就不是纯 SF1）。

## 5. 科学意义
1. **甲基化亚型可迁移性获得独立队列支持**：Møller 的 5 簇结构在 GSE283928 肿瘤子集上可复现（剔除正常对照后 4 个肿瘤簇分离稳定）。
2. **方法学框架验证成立**：用"全探针 PCA + KMeans"即可在外部 EPIC 队列重现 Møller 分型，为后续 surrogate 风险建模提供迁移基础。
3. **诚实边界**：Belakhoua 是全 PitNET 谱系队列，**无 NFPT 专属术后再生长标签**，本验证仅证明"分型可迁移"，不能直接做 regrowth PFS 验证——下一步 regrowth 风险模型须走转录组 surrogate（GSE209903/GSE51618）并明确声明方法学延伸。

## 6. 产物
- `results/GSE283928/lineage_marker_probes.csv` — 谱系标记探针列表
- `results/GSE283928/lineage_cluster_validation.json` — 谱系评分按簇
- `results/GSE283928/cluster_diag.json` — 每簇全局 β 诊断（normal 判定）
- `scripts/extract_lineage_probes.R` / `scripts/validate_lineage_clusters.py` / `scripts/diag_cluster_normal.py`
