# Nature-Style Reviewer Assessment — Manuscript V1 (NF-PitNET Methylation Taxonomy)

**评审日期**: 2026-08-06
**评估对象**: `docs/manuscript-v1.md`（ARS 第一版初稿，含全部已修复项）
**实际目标期刊**: *Neuro-Oncology* (IF 13.1)，非 Nature。下方按 nature-reviewer 五轴模拟评审；Nature 的"杰出重要性"门槛不适用于本研究，已在 Significance 轴明确标注。

---

## Review setup

- **Input scope**: 完整初稿（标题→结论→参考文献→图注→6 个强制章节→Draft Metadata）+ 配套 5 主图 / 3 主表 / S1-3 / S3 / Supp Methods 实体文件。
- **Assessment boundary**: 以稿件正文 + 结果文件（JSON/CSV/PNG）为唯一依据；图内容因本环境无法读图，仅以尺寸/图注/生成脚本间接核验，颜色与视觉细节需作者在 PDF 导出后目检。
- **Shared manuscript claim summary**: NF-PitNET 是 SF1-lineage 甲基化实体（映射 GSE283928 c0），在 4 个独立队列、2 个芯片平台可复现（c0 映射均值 84.3%），且与 Møller 最高危 k4（TPIT/PIT1→c3/c4）解离；预后风险轴在谱系内细分（k3/k4/k5 vs k1/k2 HR 2.14, p=0.072）。
- **Visible evidence base**: 联合 PCA（687,463 探针/281 样本，PC1 29.2%）、最近质心投射（4 队列）、562-DMP 探针级映射（k4→c3/c4 p=1.6×10⁻³⁰）、Neou 去标签聚类（k=3, 90% NFPT→c0）、117-idat 去标签投射（86.3%）、直接预后 Cox、562-DMP LOO 恢复（99.1%/93.7%）、CMap/Enrichr（探索性）。
- **Missing materials affecting confidence**: ① 无公开数据/代码仓库 DOI（Data Availability 仍指向本地 `scripts/`、`results/`）；② 无第二个带再生长标签的外部预后队列；③ Authors/Funding/COI 为占位符。

---

## Reviewer 1 — 技术稳健性侧重 (technical soundness)

- **Overall assessment**: 方法学扎实，循环性担忧已被主动处理；但预后主张建立在单一队列、边界显著的统计量上，且可重复性链条在"公开获取"环节断裂。
- **Who would be interested, and why**: 垂体肿瘤临床研究者、神经病理/表观遗传分类研究者——他们关心 NF-PitNET 是否真有可复现的甲基化身份、以及 k4 是否等于"侵袭性"。
- **Major strengths**:
  - 通过 Neou 去标签聚类（§3.7）与 117-idat 去标签投射（§3.8）双重阻断"质心循环"质疑，设计严谨。
  - 562-DMP 探针级映射（§3.5）把"NF-PitNET = SF1-lineage"从基因均值代理提升到探针分辨率，是本研究最强的技术贡献。
  - 归一化敏感性（Table 3 / ST3）与跨平台 HM450（§3.3）共同排除预处理/平台伪影。
- **Major concerns**:
  - **R1-M1 (statistical-rigor / claim-moderation)**: 预后主主张依赖 HR 2.14 (p=0.072) 与校正 HR 2.09 (p=0.088)，二者均 >0.05。结论段（L261）称"支持甲基化分层的生物学与预后**稳健性**（robustness）"超出证据。建议改为"提示性（suggestive）"并明确标注有限效能（n=111, 28 事件）。
  - **R1-M2 (clinical-validity / reproducibility)**: 预后验证完全依赖 Møller 117 例 + 作者提供的 k-label + 再手术终点；三个公共队列无样本级再生长标签（Limitation 1 已承认）。缺少第二个独立再生长标签队列。
  - **R1-M3 (data-resource-quality / reproducibility)**: §2.9 与 Data Availability 称"deposited in scripts/ and results/ directories of the project repository"，但全文无 Zenodo/Figshare/GitHub/DOI。*Neuro-Oncology* 要求数据可独立获取，此为投稿阻断项。
- **Technical failings to fix before the case is established**:
  - R1-M1：软化"robustness"措辞，给出效能说明。
  - R1-M3： deposit 至 Zenodo/Figshare 并引用 DOI。
- **Assessment against Nature-style criteria**: 技术完整度达 Neuro-Oncology 标准；Nature 的"杰出技术新颖性"不适用（此为验证/整合研究）。
- **Recommendation posture**: 修改后接收（minor revision 级别）；R1-M3 为硬门槛。
- **R1-m1 (writing-clarity / statistical-rigor)**: §3.2 "p approximately 0" 不精确；应给实际置换 p（如 p<1×10⁻⁴）或引用 §3.1 的精确零分布。

---

## Reviewer 2 — 原创性与科学重要性侧重 (originality / significance)

- **Overall assessment**: 概念性贡献清晰（k4 解离修正），但核心发现（NF-PitNET = SF1-lineage）与已知生物学一致，属高质量验证/整合，而非新机制或新标志物发现。
- **Who would be interested, and why**: 表观遗传肿瘤分类学读者；本研究把 Møller 2025 的分类从"单队列监督模型"提升为"多队列可独立复现的结构"，对该领域的确重要。
- **Major strengths**:
  - "k4 解离"修正是真正的智力贡献：它纠正了"最高危簇 = 侵袭性 NF-PitNET"的常见误解，且锚定 Møller 2025 的 TPIT/PIT1 定义（§4, L219-221），论证严密。
  - 跨 4 队列 / 2 平台 84.3% c0 映射 + 去标签 117-idat 86.3% 构成强证据三角。
- **Major concerns**:
  - **R2-M1 (novelty-significance)**: NF-PitNET 属 SF1-lineage 与 Sievers 2020、Capper 2018 的既有分类一致（§4, L231 已自引）。新颖性在于"首次独立多队列验证 + k4 修正"。对 *Neuro-Oncology* 足够；若投 Nature 会被 desk-reject（缺乏机制突破）。建议在引言/讨论更明确地定位为"首个独立验证"，而非暗示新发现。
  - **R2-M2 (mechanism-evidence / claim-moderation)**: §3.10 与 §4 治疗含义占显著篇幅，但证据为计算性 CMap/Enrichr（BET 已标注 LINCS 偏倚），lineage→药敏为假说。尽管有重警告，仍建议进一步压缩或明确"假说，非预测性声明"，避免权重超过证据。
  - **R1-M1 (再强调)**: 摘要（L56）"predicted regrowth (HR 2.14, p=0.072)"与结论"robustness"之间存在强度梯度，显著性边界会削弱"重要性"叙事。
- **Assessment against Nature-style criteria**: 重要性达领域级（field-significant），未达"广泛科学重要性"（broad importance）。
- **Recommendation posture**: 修改后接收；R2-M1 为 framing 建议，非硬伤。

---

## Reviewer 3 — 跨学科读者兴趣与可读性侧重 (interdisciplinary interest / readability)

- **Overall assessment**: 对专业读者清晰，Key Terms 侧栏有效；但治疗含义段落对非表观遗传读者仍偏重，且可重复性（公开仓库）缺口影响"可被复用"这一跨学科价值。
- **Who would be interested, and why**: 计算生物学/药物重定位研究者会关注 CMap 部分；但需明确这是假说。
- **Major strengths**:
  - Key Terms 侧栏（L64-70）明确定义 SF1/TPIT/PIT1、k1-k5、c0-c4、EPIC/HM450，显著降低专业门槛。
  - 图注自包含（Figure 1-5 + S3），非专家可凭图理解主逻辑。
- **Major concerns**:
  - **R3-M1 (clinical-validity)**: 再手术/再干预作为"再生长替代终点"是临床决策而非纯生物学终点，可能受术者驱动；应在局限中说明该替代性质与敏感性。
  - **R3-M2 (readability)**: 治疗含义（§3.10 + §4）对非表观遗传读者的信号噪声比偏高；建议保留假说框架但精简，突出"需 NFPA 原代细胞功能验证"这一明确下一步。
  - **R1-M3 (reproducibility → 可复用性)**: 无公开仓库使跨学科读者无法复现/复用，直接削弱"跨学科价值"主张。
- **R3-m1 (writing-clarity)**: 摘要"predicted regrowth"与结论"robustness"在强度上不完全对齐（见 R1-M1），建议统一为"提示性"口径。
- **Assessment against Nature-style criteria**: 可读性达目标期刊标准；跨学科广度有限（本质为垂体肿瘤亚专科）。
- **Recommendation posture**: 修改后接收。

---

## Cross-review synthesis

- **Consensus strengths**:
  - 562-DMP 探针级映射（§3.5）是核心技术与概念贡献。
  - 去标签聚类 + 去标签 117-idat 投射主动消解循环性质疑（experimental-design 轴强）。
  - 跨 4 队列 / 2 平台一致性 + 健康垂体内对照构成强证据三角。
- **Consensus technical risks** (≥2 位评审独立提出):
  - **C1 — 预后主张边界显著**: R1-M1（技术）+ R2（重要性叙事）+ R3（终点替代性）共同指向 HR 2.14/2.09 的 p=0.072/0.088 不足以支撑"稳健性"；需软化措辞并补效能说明。
  - **C2 — 无可复用公开仓库**: R1-M3（数据资源）+ R3（可复用性）共同指出 Data Availability 缺 DOI，为投稿硬门槛。
  - **C3 — 验证性本质限制广泛重要性**: R2-M1（新颖性）+ R3（跨学科广度）一致认为这是高质量整合/验证，非机制突破；对 Neuro-Oncology 足够，对 Nature 不足。
- **Where emphasis differs**:
  - Reviewer 1 聚焦统计边界与单一队列依赖；Reviewer 2 聚焦新颖性定位与 k4 修正的概念价值；Reviewer 3 聚焦可读性与治疗含义权重。
- **Broad-interest / significance readout**: 对垂体肿瘤表观遗传分类领域具明确价值（首个独立多队列验证 + 纠正 k4 误读）；但信号本质为亚专科，广泛科学读者兴趣有限。
- **Most important issues to resolve before a strong case is established**:
  1. R1-M3：补公开仓库 DOI（硬门槛）。
  2. R1-M1：软化"robustness"、补效能说明。
  3. R1-M2 / R3-M1：明确预后验证的单一队列依赖与再手术替代终点局限。

---

## Risk / unsupported claims

- "prognostic robustness"（结论 L261）：当前 p=0.088 不支持"稳健"，应降为"提示性/suggestive"。
- "predicted regrowth"（摘要 L56）：HR p=0.072 为边界显著，表述可接受但需与结论一致。
- CMap/Enrichr 治疗含义：全文已重警告为探索性/假说，未越界；维持当前 framing 即可。
- 无孤儿引用（1-18 全部对应），18 组关键数字均与结果文件核对一致（延续 2026-08-04 核对结论）。

---

## 到位了吗？（对用户的直接回答）

**草稿层面：已做到位。** 此前 14 项问题（2 CRITICAL + 4 MAJOR + 8 MINOR）全部修复并重新核验：
- Table 1/2 实体已生成且与正文数字完全自洽（Neou 26/32、GSE207937 34/45、Møller 101/117、GSE115783 32/34，全部对账通过）；
- 图编号重排为连续 Figure 1-5，正文 5/3/2/4/4 处引用一致，无跳号；
- `figure4_merged.png`（556 KB）已生成，Figure 4 为三面板合并图；
- §2.3 已透明披露 NN-mixing 双口径；Draft Metadata 已更新为 5,890 词 / 0 em dash；
- em dash = 0，AI 典型词 = 0。

**距投稿就绪还差 4 项快速修复（非科学硬伤）：**
1. **补公开仓库 DOI**（R1-M3，硬门槛）——Zenodo/Figshare 上传 `scripts/`+`results/`，替换 §2.9/Data Availability 的本地路径。
2. **软化"robustness"措辞**（R1-M1）——结论段鉴于 p=0.088 改为"提示性"，并补 n=111/28 事件的效能说明。
3. **替换"p approximately 0"**（R1-m1）——§3.2 给实际置换 p。
4. **填占位符**（Authors/Funding/COI）+ 明确再手术替代终点局限（R3-M1）。

这 4 项都是半天内可收口的事务性工作；科学内容本身已扎实、数字全部对齐、结构完整。
