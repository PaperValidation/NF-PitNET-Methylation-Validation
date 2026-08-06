# 经验存档：R 写出的 `.csv.gz` 实为 TSV，被 pandas 静默误读为单列

> 类别：数据 IO / 跨语言数据交换陷阱
> 首次踩坑：2026-08-01，NF-PitNET 项目 `moller_562dmp_within_lineage.py`
> 严重度：高（不报错、静默归零，极难一眼看出）

---

## 1. 症状（怎么发现它）

下游把甲基化 β 矩阵和一个探针名单做交集，结果**毫无征兆地返回 0**：

- `intersection of 0 probes`
- `0/117 samples aligned`
- 或者 `df.shape` 显示本该是 `863,306 × 117` 的宽矩阵，却变成 `863,306 × 1`。

最阴险的地方：**没有任何 Exception**，脚本照常跑完，只是所有下游结果都是空的/错的。

---

## 2. 根因

R 侧导出代码（`moller_idat_to_beta.R:87-89`）：

```r
gz <- gzfile(out_csv, "wt")                       # out_csv = "...mendeley_beta_cg.csv.gz"
write.table(cg, file = gz, sep = "\t", quote = FALSE, col.names = NA)
close(gz)
```

`write.table(..., sep = "\t")` 按**制表符**写出，但文件被命名为 `xxx.csv.gz`。
（等价陷阱：`data.table::fwrite(..., sep = "\t")` 配上 `.csv` 扩展名。）

Python 侧读取代码（初版）：

```python
reader = pd.read_csv(BETA_GZ, index_col=0, chunksize=40000)   # 默认 sep=','
```

pandas 默认逗号分隔，于是整行（行名 + 117 个样本值）被当成**一整列的字符串**读入，
索引名变成 `"cg00000001\t0.51\t0.62\t..."`。这个"假探针名"和任何真实 cg 名单都匹配不上，
交集自然为 0，且不会抛错。

> 一句话：**R 的 `write.table(sep="\t")` 导出的是 TSV，文件名却写着 `.csv`（或 `.csv.gz`），pandas 用逗号去读就会整表塌成单列。**

---

## 3. 诊断（3 秒定位）

- 看首行原始字节：`head -c 500 file | cat -A` —— 如果字段之间是 `^I`（tab）而非逗号，就是 TSV。
- pandas 只读前几行看列数：`pd.read_csv(path, nrows=3).shape[1]` —— 宽矩阵若读成 `1` 列，100% 是分隔符问题。
- 健壮嗅探（仅用于头部检查，大文件慢）：`pd.read_csv(path, sep=None, engine="python")`。

---

## 4. 根治 / 预防

**R 侧（治本）：**
- `.csv` 就用 `write.csv()`（逗号分隔）；
- 若确实要 TSV，文件名请写成 `.tsv.gz`，或写入时显式 `sep="\t"` 并在文件名体现；
- 跨语言交换宽矩阵，优先考虑 `saveRDS()` / `.parquet` / `.feather`，比纯文本更不容易踩分隔符坑。

**Python 侧（兜底）：**
- 对"宽矩阵"类文件，统一先嗅探分隔符再读；
- 解析 GEO 或 R 导出的 `.csv.gz` 时，**优先假设 TSV**：`pd.read_csv(path, sep="\t", index_col=0)`。

**修复后的读取（现场）：**

```python
# WRONG（静默零匹配）:
reader = pd.read_csv(BETA_GZ, index_col=0, chunksize=40000)
# RIGHT:
reader = pd.read_csv(BETA_GZ, sep="\t", index_col=0, chunksize=40000)
```

---

## 5. 相关脚本（已入档）

- `experience-archive/moller_562dmp_within_lineage.py` —— 本坑的发生现场与修复（提取 562 DMP、分块读 800MB TSV β 矩阵）。
- `../scripts/moller_idat_to_beta.R:87-89` —— 写出 TSV 的源头（`write.table(..., sep="\t")` → `.csv.gz`）。

---

## 6. 可复用的检查清单（每次读 R 导出的"csv"都过一遍）

1. 文件名是 `.csv` 但来源是 `write.table` / `fwrite`？→ 十有八九是 TSV。
2. 读进来先打印 `df.shape` 和 `df.columns` 前几项的原始字符串。
3. 宽矩阵列数异常（=1 或远小于预期）？→ 立刻试 `sep="\t"`。
4. 不要依赖"不报错 = 读对了"——这种陷阱恰恰不报错。
