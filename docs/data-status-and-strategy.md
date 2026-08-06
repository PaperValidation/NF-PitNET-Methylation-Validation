# NF-PitNET 数据可得性终评与策略定稿

> 更新：2026-07-29（第二轮核验后修正）

## 一、关键数据可得性终评

| 队列 | 类型 | 样本 | 公开性 | 标签/终点 | 落地状态 |
|---|---|---|---|---|---|
| **Møller 2025 (noaf269, PMID 41369112)** | EPIC 850K 甲基化 | 117 NFPitNET + 影像随访 | **原始 β 矩阵不公开**；补充材料 zip 被 reCAPTCHA 拦截（已尝试 PMC/EuropePMC 4 个端点均失败） | regrowth/PFS（SF1-lineage） | 仅全文 + 关键基因/5 簇结构 |
| **Silva-Junior 2023 (GSE207937 + GSE209903)** | EPIC 甲基化 + 配对 RNA-seq | 77（甲基化）/ 65（RNA-seq） | 公开 | subtype(GH/ACTH/NS/NA)；NS=非分泌(NFPT候选) | **已下载** matrix + samples CSV |
| **Belakhoua 2025 (noaf109, PMID 40295206, GSE283928)** | EPIC 850K 甲基化 | 118 全 PitNET 谱系 | 公开（processed matrix 874MB） | Ki-67/临床病理亚型（非 NFPT 专属 regrowth） | **已下载** 727,986 探针 × 236 列 |
| **Neou 2020 (E-MTAB-7762)** | EPIC + RNA-seq | 46 | ArrayExpress 公开 | 谱系分型 | 待下载 |
| **GSE51618** | 微阵列转录组 | 10（侵袭3/非侵袭4/正常3） | 公开 | 侵袭/非侵袭 NFPA | **已下载** |
| **GSE201439** | RNA-seq | 25 | 公开（仅 RAW.tar） | PitNET | 待下载（原始，优先级低） |

## 二、环境结论

- **WSL Ubuntu-24.04**：R 4.3.3 + Python 3.12.3 + RTX 3070 8GB，20 核，可用内存 ~4.7GB
- **代理**：Windows 主机 `127.0.0.1:8428` 通 NCBI（200）；但 WSL 跨命名空间连不到（代理只绑 Windows loopback）
- **数据流定式**：Windows 侧代理下载 → 落盘 `G:` → WSL 经 `/mnt/g`（9p 可读写）离线分析（**WSL 不需联网**）
- **下载器**：`scripts/dl_geo_v2.py`（socks5 + rdns=True + Range 头）解决 NCBI 403/rdns 问题

## 三、R 生信栈状态

已就绪：GEOquery/limma/edgeR/DESeq2/sva/WGCNA/glmnet/survival/fgsea/Seurat/SingleR/minfi/DMRcate/ggplot2/pROC/data.table/tidyverse
缺失（不影响主线，可后补或用替代实现）：clusterProfiler（用 fgsea+自写 GO/KEGG）、survminer/timeROC/cmprsk/rms（用 survival+glmnet+ggplot2 自实现 Cox/分层ROC）、ChAMP（minfi 已够）

## 四、策略定稿（修正版）

**放弃**：精确 562 DMP 探针列表（被 reCAPTCHA 硬墙拦截，与 "原始数据不公开" 一致）

**采用**（不依赖 562 精确探针，以 Møller 方法学框架 + 关键基因为锚点）：

1. **甲基化亚型可迁移性验证**（核心科学问题）
   - 在 GSE283928（Belakhoua 118 例 EPIC）上做无监督共识聚类（ConsensusClusterPlus）
   - 用 Møller 5 簇关键基因签名（CCDC141/ATP1A3/ANK2/CASC15/HDC/NCKAP5/GUCA1C/ATP2C1 等）做监督验证
   - 内存优化：Python 抽稀到关键基因 + 降维，或 R 里只保留与 EPIC v1 重叠探针

2. **转录组 surrogate 风险模型**（NF-PitNET 专属）
   - GSE209903（65 NFPT RNA-seq）+ GSE51618（侵袭/非侵袭）建 SF1-lineage 表达签名
   - 差异分析（DESeq2/limma）+ 通路富集（fgsea + 自写 GO/KEGG 映射）
   - LASSO/elastic net 特征选择 → 风险评分

3. **术后再生长主终点**：诚实声明 NFPT 专属 regrowth/PFS 公开标签稀缺（领域空白=故事点）；
   用 Møller 117 例外部验证队列的预后分离结论作为生物学支撑，在公开队列上做 surrogate 验证

4. **单细胞/空间映射**（若 Nature Commun 2025 图谱可得）：SF1 肿瘤细胞 / SPP1+ TAM / stromal niche

5. **药物重定位**：LINCS/CMap 负连接（CDK4/6、PI3K/mTOR、HDAC；SPP1-ITGAV/ITGB1 轴）

## 五、已落盘文件清单

- `G:/.../NF-PitNET/GSE51618_series_matrix.txt.gz` (1.9MB)
- `G:/.../NF-PitNET/GSE209903_series_matrix.txt.gz` (5.7KB, 仅 metadata)
- `G:/.../NF-PitNET/GSE207937_series_matrix.txt.gz` (490MB)
- `G:/.../NF-PitNET/GSE283928_MatrixProcessed_Normalized.txt.gz` (874MB)
- `G:/.../NF-PitNET/samples/*_samples.csv` (3 个队列样本元数据)
- `G:/.../NF-PitNET/metadata/master_metadata.csv` (152 样本统一表)
- `G:/.../NF-PitNET/pmc13070495_printable.html` (命门论文全文)
- `scripts/`: dl_geo.py, dl_geo_v2.py, parse_matrices.py, build_master_metadata.R, wsl_r_pkgs.R, wsl_install_*.R, inspect_gse207937.R

## 六、下一步执行顺序

1. **P1（当前）**：GSE51618 + GSE209903 转录组差异分析（小、内存无虞）
2. **P2**：GSE283928 5 簇复现（内存优化版）
3. **P3**：Neou 2020 (E-MTAB-7762) 下载 + 整合
4. **P4**：风险模型 + 外部验证
5. **P5**：单细胞/空间 + 药物重定位 + 稿件
