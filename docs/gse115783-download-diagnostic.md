# GSE115783 外部验证数据下载 — 实测诊断

> **✅ 已解决（2026-08-01 同日）**：NCBI 本环境阻断仍在，但用户改为**手动下载** `GSE115783_series_matrix.txt.gz` + `GSE115783_RAW.tar` 并放入 `E:/Medical Data/Neurosurgery/NF-PitNET/GSE115783/`。series_matrix 数据表即 ChAMP 处理 β 值，直接读入投影即可。结果：**34 NFPA 肿瘤 32/34=94.1%→c0，6 正常对照仅 33.3%→c0（特异性对照）**；跨平台 HM450 缺口已闭合，manuscript §3.9 + Figure 6 已写入。下方为当时"本环境无法自动下载"的诊断记录，保留供复盘。

> 日期：2026-08-01
> 背景：manuscript 的 Future directions (ii) 把 GSE115783（34 NFPA, Illumina 450K）列为唯一独立的 NFPA-only 甲基化阵列，用于跨平台（HM450 vs 我们的 EPIC-850K）确认 86%→c0 映射。用户指出 `Guidiance/lessons-learned/ncbi-pmc-pow` 方案可解决 NCBI 403，要求尝试。

## 结论速览

**用户给的 `ncbi-pmc-pow` 方案无法下载 GSE115783，根因是两套不同的 NCBI 子系统：**

| 维度 | `ncbi-pmc-pow` 方案 | GSE115783 所需 |
|---|---|---|
| 目标主机 | `pmc.ncbi.nlm.nih.gov/articles/instance/<pmcid>/bin/*.zip` | `ftp.ncbi.nlm.nih.gov/geo/series/.../GSE115783_series_matrix.txt.gz` |
| 反爬机制 | **SHA256 前缀 PoW 挑战页** + 种 cookie 放行 | **IP 级封禁**（数据中心 IP 的 GET 一律 403/超时） |
| 适用场景 | PMC 补充材料二进制 | GEO 系列 β 矩阵 |

二者主机、机制、解法完全不同。PoW 解法对 GEO FTP 无效。

> 注：用户给的目录与本项目 `docs/ncbi-pmc-pow-extraction-guide.md` **同源**（PMC SHA256-PoW 实战记录），方案本身已收录。

## 已实测的所有通道（全部失败）

| # | 通道 | 结果 |
|---|---|---|
| 1 | `ftp.ncbi.nlm.nih.gov/.../GSE115783_series_matrix.txt.gz` 直连 | **HTTP 403**（HEAD 可过、GET 被拒）/ 偶发超时 — 数据中心 IP 被封 |
| 2 | `www.ncbi.nlm.nih.gov/geo/download/?acc=...&file=matrix/...` | **HTTP 404** — 端点/相对路径不对 |
| 3 | `www .../acc.cgi?acc=GSE115783&form=text&view=full` | 仅返回 **3894B 元数据 SOFT**，**不含 β 矩阵**（矩阵只在独立 matrix 文件里） |
| 4 | HTTP 代理 `127.0.0.1:8428` → ftp/矩阵 | NCBI **直接关连接**（代理出口 IP 也被封） |
| 5 | SOCKS5 `127.0.0.1:1987` | **未监听**（WinError 10061）— 本会话代理没开 |
| 6 | E-utilities `eutils.ncbi.nlm.nih.gov` | **可达**（拿到 GDS UID=200115783），但 efetch 不返回原始 β 矩阵 |
| 7 | EBI ArrayExpress `www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-115783/` | **页面 200 可达**（EBI 确实镜像了 GSE115783） |
| 8 | EBI 文件直链 / BioStudies API `api/v1/studies/E-GEOD-115783/files` | **404 / API 404** — 矩阵经 JS/BioStudies 前端暴露，简单直链取不到 |
| 9 | EBI FTP `ftp.ebi.ac.uk/.../E-GEOD-115783/` | **404** — 归档路径已迁移 |

## 唯一真实线索 & 可行出路

- **EBI / BioStudies 确实托管 E-GEOD-115783**（步骤 7 实测 200），但数据文件通过 BioStudies 的 JS 前端 + API 暴露，脚本化直链拿不到。
- **可行出路（任选其一）**：
  1. **经 EBI 网页 UI 手动下载**：浏览器打开 `https://www.ebi.ac.uk/biostudies/studies/E-GEOD-115783`，手动导出 β 矩阵/处理数据，丢到 `E:/Medical Data/Neurosurgery/NF-PitNET/GSE115783/`，我再跑投影。
  2. **等 NCBI 访问恢复**：换住宅 IP 的机器、或本机代理（8428/1987）重启后重试 `ftp.ncbi.nlm.nih.gov` 直连。
  3. **换一个本环境可达的 HM450 NFPA 队列**（若接受替代而非 GSE115783 本身）。

## 对 manuscript 的影响（无伤）

- 三队列 **EPIC-850K** 独立外部验证（Neou 81% / GSE207937 76% / Møller 86%，均值 81%）**不受影响、仍然成立**。
- 跨平台（HM450）维度缺口已如实写入 Limitations / Future directions (ii)，标注为"designated next step, currently blocked"。
- 不夸大、不编造 GSE115783 结果。

## 探测脚本（留存备查）

- `scripts/probe_ncbi_endpoints.py` / `probe_ncbi_endpoints2.py` / `probe_ncbi_endpoints3.py` — 各端点吞吐探测
- `scripts/probe_ncbi_proxy.py` — 直连 vs 8428 代理对比
- `scripts/probe_ncbi_1987_softfull.py` — 1987 SOCKS5 + SOFT full 长超时尝试
- `scripts/find_gse115783_matrix_url.py` — 抠 GEO 页面真实矩阵链接
- `scripts/update_memory_ncbi.py` — 本次结论写入记忆
