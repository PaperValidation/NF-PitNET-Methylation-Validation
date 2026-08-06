# 叙事重构方案：基于 cluster-风险映射结果

> 对应预审委员会致死级问题 2
> 日期：2026-07-31

> ⚠️ **已被 562-DMP 探针级验证 supersede（2026-07-31 晚）**：本文基于 `cluster_risk_mapping.py` 得出的核心结论 **"GSE283928 c0 ≈ Møller k4（HR 3.63，最高风险）" 已被证伪**。
> 决定性证据见 `docs/moller-dmp-probe-validation.md`：562 个精确 DMP 的探针级方向一致性映射显示 **c0 对应 Møller 的 SF1-lineage 簇（k1/k2/k3/k5）**，而最高风险 **k4（TPIT/PIT1, HR 3.63）映射到 c3/c4（p=1.6e-30）**。下文"核心发现"与"方案 A"中"c0 = k4 / 最高风险"的表述均已失效，请以探针级验证报告为准；manuscript 已同步订正。

---

## 核心发现：预审质疑的答案是正向的

预审委员会提出的最致命质疑是："NFPT 映射到 GSE283928 cluster 0，但 cluster 0 对应 Møller 的高风险簇还是低风险簇？如果是低风险，核心叙事方向就反了。"

**审计结果（`cluster_risk_mapping.py`）：**

| GSE283928 簇 | N | Møller 关键基因 β | vs c0 perm-p | Lineage Score | 风险对应 |
|---|---|---|---|---|---|
| **c0** | 70 | **0.6014** (SF1-lineage 类比最高) | — | -0.1073 (SF1-lineage 富集) | **SF1-lineage 类比（k1/k2/k3/k5）；k4→c3/c4（见 562-DMP 探针级验证，已 supersede 本表原"≈k4"结论）** |
| c2 | 2 | 0.6175 | ns (n=2) | +0.0234 | 不可靠 |
| c3 | 22 | 0.5834 | **0.0009** | +0.0163 | ≈ k3/k5（SF1, HR 2.57/2.69） |
| c4 | 24 | 0.5776 | **0.0001** | -0.0581 | ≈ k1/k2（SF1, 低风险参考） |

**三个关键结论（基于 `cluster_risk_mapping.py`，⚠️ 已被 562-DMP 探针级验证修正）：**

1. **Cluster 0 在 3 个主要簇中 Møller 关键基因甲基化最高**（c0 > c3 > c4, p < 0.001）
   - 但此信号反映 **SF1-lineage 谱系身份**（Møller-key 探针集含 *NR5A1*/SF1），而非 k4 风险
   - → 探针级验证表明 c0 对应 Møller 的 **SF1-lineage 簇（k1/k2/k3/k5）**

2. **Cluster 0 的 lineage score 最低（-0.1073）**——SF1-lineage 富集，*而非* TPIT/PIT1-like
   - ~~Møller 的 k4 是唯一的 TPIT/PIT1 簇，HR 3.63（最高风险）→ c0 ≈ k4~~ **（错误，已证伪）**
   - 探针级映射显示 **k4（HR 3.63, TPIT/PIT1）映射到 c3/c4**（p=1.6e-30），与 c0 完全反向

3. **NFPT→c0 富集是真实的**（81.2% vs null 32.5%, perm p ≈ 0, LOO accuracy 89.8%）
   - 不是"大簇引力场"伪影；NFPT 确实稳定映射到 SF1-lineage 类比簇

**→ 修正后的叙事：NFPT 映射到可跨三队列复现的、独特的 SF1-lineage 甲基化身份（c0 = Møller k1/k2/k3/k5 类比）；Møller 最高风险 k4（TPIT/PIT1）对应分泌型组（c3/c4），而非 c0。原"NFPT → 最高风险 k4"卖点被证伪，但换得更稳健的"分型框架探针级可复现"结论。详见 `docs/moller-dmp-probe-validation.md`。**

---

## 方案 A（推荐）：强化当前叙事

**核心论点升级**：

原叙事："NFPT 映射到 GSE283928 cluster 0"（模糊——cluster 0 的风险等级未知）

升级为："NFPT 跨三个独立 EPIC 队列稳定映射到甲基化最高风险簇——该簇具有最高的 Møller 关键基因甲基化（p < 0.001 vs 其他主要簇）和最 TPIT/PIT1-like 的谱系特征，对应 Møller 2025 的 k4（HR 3.63，术后再生长最高风险）。"

**新增证据链条（本次审计产出）：**

```
NFPT (3 队列一致) → GSE283928 c0
                    ↓ (cluster_risk_mapping.py)
                    c0 = 最高 Møller key β (0.6014, p<0.001 vs c3/c4)
                    c0 = 最低 lineage score (-0.1073, TPIT/PIT1-like)
                    ↓
                    ≈ Møller k4 (HR 3.63, 最高再生长风险)
                    ↓
  NFPT → 最高风险甲基化亚型 → 术后再生长风险最高的生物学亚群
```

**叙事强化要点：**
1. 不再只说"可迁移性"，而说"NFPT 稳定映射到最高风险亚型"
2. 引入 lineage 维度：NFPT 中非 SF1 主导的亚群可能携带最高再生长风险
3. 与 Møller 的 k4（HR 3.63, 唯一非 SF1 簇）形成直接对应

**适合期刊**：Neuro-Oncology（IF 13.1）可尝试；Journal of Neuro-Oncology 稳妥

**风险**：GSE283928 是混合谱系队列，c0 的 lineage score 最低可能反映样本组成（c0 可能混入了更多 TPIT/PIT1 谱系肿瘤），而非 NFPT 本身的谱系特征。需要在 Discussion 中讨论这一混杂。

---

## 方案 B（保守）：降维为分类学可重复性研究

如果方案 A 的 lineage 混杂问题无法解决，退而求其次：

**核心论点**："在不使用 Møller 分类器参数的情况下，完全从零开始在三独立 EPIC 队列上复现了 NF-PitNET 甲基化亚型结构，证明该分型是数据驱动的普适生物学现象。"

**优势**：
- 不需要建立 cluster ↔ risk 的对应关系
- 不需要 regrowth 标签
- 纯方法学验证，审稿人不会质疑"没有预后终点"

**劣势**：
- 增量贡献偏薄（Møller 已做了外部验证）
- 丢失了"最高风险"这个叙事亮点
- 更适合 Journal of Neuro-Oncology 而非 Neuro-Oncology

---

## 方案 C（最强但需额外工作）：加入 Neou 内部亚型分析

在方案 A 基础上，对 Neou 86 个 EPIC 样本做**独立无监督聚类**（k=5），看是否能发现与 Møller 5 簇对应的亚结构。如果能：
- 直接在 Neou 内部建立 cluster ↔ Møller 簇的对应
- 不再依赖 GSE283928 的 cluster 作为中间桥梁
- 增量贡献最大：第三个独立队列的完全从头分型

**所需工作**：~2 小时（Neou β 矩阵已有，需跑 ConsensusClusterPlus 或 KMeans + 注释）

---

## 推荐：方案 A + C 组合

1. **立即执行方案 A**：用已有的 cluster-risk 映射结果强化叙事
2. **追加方案 C**：对 Neou 做独立聚类，作为 supplementary analysis
3. **在 Discussion 中诚实讨论**：
   - c0 的 lineage 混杂（GSE283928 是混合谱系队列）
   - Møller key 评分为基因级均值，非 DMP 探针级
   - 无样本级 regrowth 标签 → 风险对应为生物学推断

这样既不丢失叙事力度，又保持了方法学诚实。
