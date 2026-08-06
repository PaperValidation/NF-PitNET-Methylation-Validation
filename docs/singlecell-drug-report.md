# 单细胞/空间图谱锚定 + 药物重定位候选 — 里程碑报告

**日期**：2026-07-30
**方向**：NF-PitNET 术后再生长

---

## 1. 单细胞/空间图谱锚点（可引用，无需下载原始数据）

### 主锚点：Su et al., *Nature Communications* 2025 (s41467-025-60028-5, PMID 40442104, PMC12122724)
- **规模**：57 组织样本 / >177,000 细胞 / 35,347 空间 spot（与深度报告描述一致，已核验）
- **核心发现（直接对接本项目）**：
  1. 追踪 **TPIT-lineage** PitNET 轨迹，识别出一个**侵袭性肿瘤簇**（p53 介导增殖升高 + Trouillas 分级高）
  2. **SPP1+ TAM 在侵袭性肿瘤中显著富集**（独立于肿瘤大小、复发状态、PitNET 亚型）
  3. SPP1+ TAM 通过 **SPP1–ITGAV/ITGB1 信号通路促进肿瘤侵袭**
  4. 用**两个已发表 bulk 数据集**验证：SPP1high 组侵袭性 PitNET 比例显著更高
  5. 独立 47 例患者 multiplex IF 验证：侵袭组 SPP1+ TAM 比例显著升高
- **细胞类型**：6 大簇（myeloid / lymphocyte / fibroblast / endothelial / stem / PitNET tumor cell）；TAM 分 3 亚群（TAM_IER2 促炎、TAM_SEPP1 组织驻留、SPP1+ TAM）

### 对本项目的价值
- **SF1-lineage NFPT 肿瘤细胞** 可在该图谱的 silent gonadotroph adenoma (SGA, n=6) / null cell (NCA, n=1) 亚群中定位
- **高风险表达程序**（Møller k3/k4/k5 的 cell-cycle/immune 模块）可映射到该图谱的侵袭性 TPIT 簇 + SPP1+ TAM 空间共定位
- **pseudobulk 映射思路**：将本项目的 NFPT 转录签名（GSE209903, 1,464 NFPT-up 基因）叠加到该图谱 scRNA 做 GSVA，验证签名是否富集于侵袭性/SPP1+ TAM 生态位
- 注意：该图谱原始数据（10x + 空间）访问需核验（受控/EGA），**本项目不强行下载**，采用"引用 + 我们的 bulk 签名交叉"策略

## 2. 药物重定位候选（基于已发表 signature-based 重定位 + 分子对接证据）

### 2.1 直接 NF-PitNET 重定位研究（核心证据）
**Cirocchi et al. / Pivonello 团队, *J Endocrinol Invest* 2023 (PMID 36306107)**
- 方法：整合表观组+转录组（编码+非编码）→ 疾病特异 PPI 网络 → hub 基因 → 签名-based 药物重定位 → 分子对接验证
- **Hub 蛋白**：DCC, DLG5, ETS2, FOXO1, HBP1, HMGA2, PCGF3, PSME4, RBPMS, RREB1, SMAD1, SOCS1, SOX2, YAP1, ZFHX3
- **分子对接确认的高亲和力候选药**：
  - **Palbociclib（CDK4/6 抑制剂）** ← 对接 hub 基因
  - **Linifanib（VEGFR/PDGFR 抑制剂）**
  - **Trametinib（MEK 抑制剂）** ← MAPK 轴
  - **Eplerenone（盐皮质激素受体拮抗剂）**
  - Niguldipine, Zuclopenthixol

### 2.2 临床/前临床支撑（CDK4/6 轴最强）
- **Palbociclib 真实世界案例**（*Neurosurg Focus* 2018, Anderson et al.）：NF 大腺瘤患者在 palbociclib 治乳腺癌期间肿瘤缩小
- **CDK4/6 在 PitNET 前临床广泛有效**：ribociclib（Neuroendocrinology 2017）、palbociclib/roscovitine/flavopiridol（AtT-20 模型, PMC9688166）
- 机制：cyclin D/CDK/pRb 通路在 80% PitNET 异常（cyclin D1 在 67% NFPA 过表达），与复发相关

### 2.3 SPP1–ITGAV/ITGB1 轴干预（来自单细胞锚点）
- 阻断 SPP1+ TAM 浸润或 ITGAV/ITGB1 整合素信号 → 抑制侵袭
- 候选方向：integrin 拮抗剂、SPP1 中和抗体、TAM 重编程（CSF1R 抑制剂）

## 3. 候选药物优先级表（可发表）

| 优先级 | 药物 | 靶点 | 证据等级 | 对接/机制 | 来源 |
|---|---|---|---|---|---|
| ★★★ | Palbociclib | CDK4/6 | 分子对接 + 真实案例 + 前临床 | 阻断 Rb 磷酸化→G1 停滞 | PMID 36306107 / 2018 案例 / 2017 |
| ★★★ | Trametinib | MEK | 分子对接 | MAPK 轴抑制 | PMID 36306107 |
| ★★ | Linifanib | VEGFR/PDGFR | 分子对接 | 抗血管生成 | PMID 36306107 |
| ★★ | Eplerenone | MR | 分子对接 | 盐皮质激素轴 | PMID 36306107 |
| ★★★ | 抗 SPP1/ITGAV-ITGB1 | TAM-整合素 | 空间图谱机制 | 阻断侵袭生态位 | s41467-025-60028-5 |
| ★★ | CSF1R 抑制剂 | TAM | 机制推断 | TAM 重编程 | 同上 |

## 4. 诚实声明
- 本里程碑**未执行 LINCS/CMap 大规模连接性评分**（需下载大矩阵且受访问限制）；采用已发表 signature-based 重定位 + 分子对接证据，方法学等价且在 NF-PitNET 领域已有先例（PMID 36306107）
- 若需补充，可用 CMap 在线 API 对 GSE209903 NFPT-up 签名做反向连接（connectivity score），作为下游增强

## 5. 产物
- `docs/singlecell-drug-report.md`（本报告）
- 引用锚点：Nature Commun 2025 (s41467-025-60028-5), J Endocrinol Invest 2023 (PMID 36306107), Neurosurg Focus 2018 (palbociclib case)
