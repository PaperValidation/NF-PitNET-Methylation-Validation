# [AUTHOR_ACTION] 待填清单 — 需要你提供的信息

> 稿件 `docs/manuscript-v2.md` 中所有需要临床/作者团队提供真实信息的占位。除这些外，稿件已是 Neuro-Oncology 投稿格式终稿（16 条文献已核验、格式模块全部达标）。

## 1. 作者与通讯信息（稿件 第 18、22 行）

- [ ] **作者列表**（line 18）：姓名、单位、ORCID（Neuro-Oncology 要求每位作者 ORCID）
- [ ] **通讯作者**（line 22）：姓名 + 完整联系（单位、地址、邮箱、电话）+ ORCID
- [ ] 确认作者顺序与贡献匹配（与 CRediT 一致）

## 2. CRediT 贡献声明（稿件 第 360–364 行）

- [ ] 把 `[Author 1]` / `[Author 2]` / `[Corresponding author]` 占位替换为真实姓名
- [ ] 按实际贡献勾选角色（建议骨架已写好：分析/写作 vs 监督/资源）

## 3. 利益冲突（line 368）

- [ ] 全体作者确认 "no competing financial or non-financial interests"（如有冲突须逐条列出）

## 4. 基金（line 372）

- [ ] 填写实际资助（基金号 + 资助机构）；若无资助，用模板句：
  > "This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors."

## 5. 数据/代码归档联动（随作者名单一起改）

- [ ] **Zenodo**（DOI 10.5281/zenodo.21840024）：把 creators 占位 "Tidus to be replaced" 替换为真实作者名单（与稿件一致）；接收后点 **Publish**
- [x] **GitHub**（github.com/PaperValidation/NF-PitNET-Methylation-Validation）：**已是 public**（2026-10-18 核实），投稿可直接引用该链接
- [ ] 稿件 Data Availability（line 352）里的 Zenodo DOI 已就位，无需再改（仅确认文字）

## 6. 安全事项（与投稿无关，但必须处理）

- [ ] **撤销 GitHub token**：github.com/settings/tokens（`ghp_o42hT…` 旧 token；**另有 `ghp_ecUct…` 本次推送用，两个都要撤**，用完即撤）
- [ ] **撤销 Zenodo token**：zenodo.org/account/settings/applications/tokens/（`i7xQtIjp…`）

## 7. 附：当前可勾选已完成项

- ☑ 文献 16 条全部真实 + 被正文引用 + 编号连续（nature-ref-verifier 复核，commit cf27614）
- ☑ Neuro-Oncology 格式模块全部达标（摘要 226 词 / Keywords 5 / Key Points 3 / Importance 114 词 / Running title 45 字符）
- ☑ Cover letter 草稿已写：`docs/cover-letter.md`（只差填第 1、2 项信息）
- ☑ 投稿前终检清单已写：`docs/submission-checklist.md`
- ☑ 数据路径规范化完成（commit b5d7fd5 + 后续；脚本经 `_paths.py`/`paths.R` 读 .env）
- ⏳ 待办：作者/通讯/基金/COI（即本清单 1–4 项）+ token 撤销（第 6 项）
