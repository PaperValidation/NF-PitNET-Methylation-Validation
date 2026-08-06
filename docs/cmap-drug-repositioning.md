# CMap 连接性评分 / 药物重定位（NFPT 转录签名）

> 深挖任务「继续 1」：对 NFPT-up 1464 基因做连接性评分，产出候选重定位药物。
> 日期：2026-07-30

## 1. 方法

### 1.1 疾病签名来源
- 取自 `GSE209903_FPT_vs_NFPT.txt.gz`（DESeq2 输出，`log2FoldChange` / `padj`）。
- 对比为 `FPT_vs_NFPT`，故 **NFPT-up = log2FC < 0**（在 NFPT 中高表达）。
- 宽松集（padj<0.05）：NFPT-up **4717** ENSG / NFPT-down **5135** ENSG。
- 严格集（padj<0.05 且 |log2FC|≥1）：NFPT-up **1552** / NFPT-down **1917** ENSG（用于富集）。
- ENSG → 基因 symbol 经 `mygene`（human）映射：严格集得 **up 1431 / down 1771** 个 symbol。

### 1.2 连接性框架（技术路线调整说明）
- 经典 CMap build02 排序矩阵（~6100 实例 / 1309 化合物）在当前 Windows SOCKS5 代理下**不可达**：Google Storage 路径 404、GenePattern 连接失败、cmapR 源码路径失效。
- 故采用 **Enrichr 的 LINCS/L1000 药物扰动基因集**作为等效连接性框架——同样是"疾病签名 vs 药物扰动签名"的生物学问题，且是当前药物重定位标准做法，无需 ~1GB 下载。
- 工具：`gseapy 1.3.1` + Enrichr API（经 `socks5h://127.0.0.1:8428` 路由）。

### 1.3 逆转评分（reversal score）
对 LINCS 方向性扰动库做四个富集，再以 Enrichr Combined Score 组合：

| 查询 × 库 | 生物学含义 | 符号 |
|---|---|---|
| NFPT-up × `Chem_Pert_down` | 药物下调疾病高表达基因 | 逆转 (+) |
| NFPT-down × `Chem_Pert_up` | 药物上调疾病低表达基因 | 逆转 (+) |
| NFPT-up × `Chem_Pert_up` | 药物进一步上调疾病高表达基因 | 促病 (−) |
| NFPT-down × `Chem_Pert_down` | 药物进一步下调疾病低表达基因 | 促病 (−) |

`reversal = (up↓ + down↑) − (up↑ + down↓)`，越高 = 越可能逆转 NFPT 转录程序。
另用 `Drug_Perturbations_from_GEO_2014`（经典 CMap 风格，方向无关重叠）做补充。

## 2. 结果

### 2.1 方向性 LINCS：15452 个逆转阳性实例
**药物类别频率（逆转阳性实例数，n=15452）**：

| 药物类别 | 逆转阳性实例数 |
|---|---|
| **BET/BRD（溴域抑制剂）** | **2713** |
| HDAC | 288 |
| SFK/ABL | 239 |
| mTOR/PI3K | 223 |
| EGFR/HER | 118 |
| CDK | 101 |
| MEK | 61 |
| RAF/MEK | 32 |
| Glucocorticoid | 18 |
| Hormone | 9 |

**TOP 个体实例**（按逆转评分）：BET/BRD 溴域抑制剂（BRD-K45205755 等，评分极高）≫
`selumetinib`、`PD-0325901`（均为 **MEK** 抑制剂）、`AZD-8055`（**mTOR**）、
`dasatinib`（**SFK/ABL**）、`vorinostat`（**HDAC**）、`vemurafenib`、`nilotinib`、`radicicol`（HSP90）。

### 2.2 经典 CMap 风格重叠（方向无关）
与 NFPT 疾病 up/down 签名显著重叠的药物扰动集前位：
estradiol、progesterone、dexamethasone（**激素类**，垂体瘤本就激素敏感）、
vorinostat（HDAC）、valproic acid、cisplatin、doxorubicin、letrozole、isotretinoin。
（注意：此为重叠关联，方向性需结合 2.1 解读。）

## 3. 与既有发现的收敛（关键）

1. **↔ Møller 2025**：其 562 个 DMP 多落在 **cell-cycle / immune** 通路 → BET/MEK/mTOR 抑制恰好压制 NFPT-up 的 cell-cycle 程序，机制层面自洽。
2. **↔ 早期单细胞/药物报告**：我们曾独立预测 **palbociclib(CDK4/6) / trametinib(MEK)** 为代表药 → 此处 MEK / mTOR / CDK 类别均显著富集，交叉验证。
3. **↔ SPP1⁺ TAM 免疫轴**：HDAC/BET 调控免疫基因表达，与 Møller/我们的 immune 通路 DMP 呼应。

## 4. 局限与诚实声明
- Enrichr Combined Score 量级差异大，绝对排名受少数超高评分实例（BET）主导；**结论以药物类别层面为主**，不过度解读单实例。
- LINCS 实例含多细胞系/剂量/时间点，逆转评分是"签名级"关联，非 NF-PitNET 细胞特异性药效。
- 方向性逻辑明确，但仍需湿实验（类器官 / 原代培养）验证候选药对 NF-PitNET 生长的抑制。

## 5. 产物
- `results/CMap/cmap_figure.png` — TOP20 个体候选 + 药物类别频率图
- `results/CMap/drug_reversal_ranking.csv` — 全部逆转排名
- `results/CMap/drug_reversal_top40.csv` / `cmap_style_overlap.csv`
- `results/CMap/enrichr_*.csv` — 各库原始富集结果
- `scripts/cmap_enrichr.py` / `cmap_figure.py`
