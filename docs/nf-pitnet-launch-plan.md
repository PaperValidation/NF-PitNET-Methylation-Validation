# NF-PitNET 术后再生长 — 资料摸底与启动方案

> 整理日期：2026-07-29
> 范围：仅参考深度调研报告中「无功能垂体神经内分泌肿瘤（NF-PitNET）术后再生长」方向
> 数据根目录：`E:\Medical Data\Neurosurgery\NF-PitNET`（2026-10-18 更新；原 `G:\CMCC Sync\...` 已失效）
> 参考手册：`E:\Medical Data\Guidiance\{生信数据库汇总,生信期刊汇总,生信研究思路汇总}.md`

---

## 一、资料盘点（已读）

| 来源 | 关键内容 |
|---|---|
| `docs/deep-research-report.md`（第 157–197 行） | NF-PitNET 方向的临床背景、可用公开数据集表、推荐题目与路线、风险/药物重定位/期刊 |
| `生信数据库汇总.md` | 35 个库（GEO/TCGA/ChAMP/CellPhoneDB 等）访问方式、坐标版本、引用规范 |
| `生信期刊汇总.md` | 10 本纯生信 SCI 的 IF/分区/OA/APC（本方向最匹配 Neuro-Oncology 13.1、Journal of Neuro-Oncology） |
| `生信研究思路汇总.md` | 27 条套路（#4 甲基化驱动亚型、#23 ML 风险模型、#27 多模态整合直接可用） |

---

## 二、命门结论：原始数据可得性（已逐项核验）

### 2.1 命门论文（甲基化主轴）
- **Møller et al., Neuro-Oncology 2025, PMID 41369112 / DOI 10.1093/neuonc/noaf269**
- 117 例回顾性 NFPitNET + Illumina EPIC 850K + 长期影像随访；5 个甲基化簇（k1–k5：4 个 SF1 主导 + 1 个 TPIT/PIT1），**562 个 DMP**，k3/k4/k5 高风险，分类器 ~97% 准确率，3 个外部队列验证。
- **关键限制：117 例原始 β 矩阵 NOT public**（Data Availability Statement 仅声明作者持有 + 1MB 补充材料）。→ 不能直接"复现 117 例甲基化模型"。

### 2.2 但拿到了无价的锚点 + 3 个公开外部队列
- **562 DMP 官方基因集（锚点）**：从 PMC 全文（已抓至 `pmc13070495_printable.html`）提取的关键基因为：
  `CCDC141, ATP1A3, ANK2, CASC15, HDC, NCKAP5, GUCA1C, ATP2C1, RP11-497G19.3, RSPH6A, RP11-204E9.1, RP11-260M19.2, RP11-219G17.4`
  （这些可作为 surrogate 转录组建模的官方签名种子 + 富集种子）
- **三个外部验证队列（均已公开）**：
  | 队列 | 登录号 | 规模 | 备注 |
  |---|---|---|---|
  | Silva-Junior 2023 (JCEM) | **GSE207937**（甲基化，77 例含 46 NFPT）+ **GSE209903**（配对 RNA-seq） | 46 NFPT | 开放 |
  | Neou 2020 (Cancer Cell) | **E-MTAB-7762**（ArrayExpress，配对甲基化+转录组） | 46 | 开放，原始矩阵可得 |
  | Belakhoua 2025 (Neuro-Oncology noaf109) | **GSE283928** | 118 | 已公开处理矩阵（874MB），含 12 例复发标签 |

### 2.3 开放 RNA-seq（发现/验证集，来自报告）
- **GSE51618**：10 例（非侵袭/侵袭 NFPA + 正常垂体微阵列）
- **GSE201439 / PRJNA831798**：25 例 PitNET RNA-seq（含促性腺轴非功能亚型）

### 2.4 单细胞 / 空间
- 2025 Nature Communications：57 组织 / >177,000 细胞 / 35,000 spots（原始访问需核验）
- **HRA003483**：人垂体肿瘤 scRNA-seq，**受控访问**

---

## 三、比报告更优的启动策略（已修正）

报告原话：117 例拿不到就"降级为转录组代理"。实测后我们有**比降级更强的打法**：

> **主轴 = 甲基化亚型「可迁移性」+ 转录组代理风险模型 + 外部独立验证**
> 1. 以 noaf269 的 **562 DMP / 5 簇签名**为官方锚点（不是凭空建模）；
> 2. 在 **3 个公开外部甲基化队列（240+ 例）** 上复现 5 簇结构并验证预后分离（Belakhoua 含 12 例真实复发标签）；
> 3. 用已发表高风险甲基化相关基因集，在 **GSE51618 / GSE201439 RNA-seq** 上构建 **transcriptomic surrogate risk score**（诚实声明方法学延伸）；
> 4. 把高风险表达程序映射到 **2025 Nat Commun 单细胞/空间图谱**（SF1 肿瘤细胞 / SPP1+ TAM / stromal niche）；
> 5. 主终点 **postsurgical regrowth / PFS**，做 KM + 时间依赖 ROC + 校准 + DCA，独立队列验证；
> 6. 药物重定位围绕 cell cycle / adhesion / immune / SPP1-ITGAV-ITGB1 轴做 LINCS/CMap 负连接（CDK4/6、PI3K/mTOR、HDAC）。

---

## 四、环境状态（已探）

- **WSL Ubuntu-24.04 在线**；R 4.3.3、Python 3.12.3、**RTX 3070 8GB GPU** 就绪（单细胞/空间/深度学习算力够）。
- **代理**：Windows 主机侧 `socks5://127.0.0.1:8428` 访问 NCBI 返回 200（可用）；**但 WSL 跨命名空间连不到**（WSL `eth4` 在 192.168.30.x，网关 192.168.30.1=Windows 真实 IP，而代理只绑 Windows loopback）。
  - **解决**：让代理软件开启"允许局域网连接"并暴露到 `192.168.30.1:8428`；或在 WSL 内设置 `http_proxy/https_proxy=socks5://192.168.30.1:8428`。下载若仍不便，可从 Windows 侧（Git Bash + 127.0.0.1:8428）下载后落盘到 `G:\...NF-PitNET`，WSL 直接读 G: 挂载。

---

## 五、阶段化工作计划（12–14 周，可按需压缩）

| 阶段 | 任务 | 关键产出 |
|---|---|---|
| **P0 环境** | WSL 装 R 生信栈（GEOquery/limma/edgeR/DESeq2/sva/WGCNA/glmnet/survival/clusterProfiler/fgsea/Seurat/SingleR/minfi/ChAMP/DMRcate）+ Python（scanpy/scvi-tools/anndata/lifelines/gseapy/PyDESeq2）；打通代理 | 可复现环境（renv/Docker） |
| **P1 数据** | 用 GEOquery/ArrayExpress 下载 GSE51618、GSE201439、GSE207937、E-MTAB-7762、GSE283928；重建样本表（分组/平台/随访状态） | 样本元数据表 + 原始矩阵 |
| **P2 差异** | GSE51618/GSE201439 标准化、离群、limma/DESeq2 差异；clusterProfiler/fgsea 富集（cell cycle/adhesion/immune/ECM） | DEG 表 + 通路 |
| **P3 甲基化/代理** | 外部队列复现 5 簇 + 562 DMP 预后分离；RNA-seq 构建 surrogate 风险评分 | 亚型可迁移性结果 + surrogate 模型 |
| **P4 单细胞** | Nat Commun 图谱注释（SF1 肿瘤细胞/SPP1+ TAM/stromal）+ pseudobulk 映射高风险程序 | 空间微环境锚定 |
| **P5 风险模型** | 主终点 regrowth/PFS；LASSO/elastic net；KM/ROC/校准/DCA；独立队列验证 | 预后模型 + 验证 |
| **P6 转化/稿件** | LINCS/CMap 药物重定位（四步过滤）；主图/补图；初稿 | 投稿稿（Neuro-Oncology / J Neuro-Oncol） |

---

## 六、目标期刊（已核验）
- **Neuro-Oncology**：IF 13.1，混合刊（可选非 OA 无 APC），上限最高 ✅
- **Journal of Neuro-Oncology**：Springer 订阅体系，可选非 OA ✅
- 备选：Acta Neuropathologica Communications / npj Precision Oncology（多为 OA，次选）

---

## 七、风险与缓解
1. **117 例原始矩阵不公开** → 用 3 公开外部队列 + 562 DMP 锚点做可迁移性，诚实声明 surrogate 方法学延伸（审稿人接受）。
2. **单细胞原始访问需核验 / HRA003483 受控** → 优先用 2025 Nat Commun 已发表图谱；HRA003483 走受控申请或仅作参考。
3. **标签口径异质**（progression/regrowth/invasion 是否等义）→ 元数据整理置顶，明确每队列定义。
4. **WSL 代理跨命名空间** → 见第四节，必要时 Windows 侧下载落盘。

---

## 八、立即启动项（建议第一步）
- [ ] P0：打通 WSL 代理 + 安装 R/Python 生信栈（先装轻量部分验证代理可用）
- [ ] P1：优先下载体量小、开放、且带标签的 **GSE51618 + GSE283928（含 12 例复发）**，建立样本表
- [ ] 同时向 E-MTAB-7762（Neou）与 GSE207937（Silva-Junior）发起下载
