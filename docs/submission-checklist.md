# Neuro-Oncology — Pre-Submission Checklist

> Target journal: *Neuro-Oncology* (Oxford University Press, original article).
> Manuscript: `docs/manuscript-v2.md` (16 verified Vancouver references; Figures 1–5 + S3; Tables 1–3; Supplementary Tables 1–3; Supplementary Methods).

## A. Front-matter modules (all programmatically verified — re-check after any edit)

| # | Item | Spec | Current | Done |
|---|------|------|---------|------|
| A1 | Structured abstract | ≤250 words; Background ends with aim sentence | 226 words | ☐ |
| A2 | Keywords | ≤5 | 5 | ☐ |
| A3 | Key Points | 2–3 bullets, each ≤85 characters | 3 (79/71/77) | ☐ |
| A4 | Importance of the Study | ≤150 words | 114 | ☐ |
| A5 | Running title | ≤50 characters | 45 ("NF-PitNET methylation taxonomy across cohorts") | ☐ |
| A6 | Title page | Authors + affiliations + ORCIDs + corresponding author contact + word count | [AUTHOR_ACTION] | ☐ |
| A7 | Total word count | ~7,000 stated | stated | ☐ |

## B. Body & mandatory sections (all present in v2 — confirm content final)

| # | Item | Done |
|---|------|------|
| B1 | IMRaD structure (Introduction/Methods/Results/Discussion) | ☐ |
| B2 | Data Availability — GEO accessions + ArrayExpress + Mendeley + Zenodo DOI 10.5281/zenodo.21840024 | ☐ |
| B3 | Ethics Declaration — secondary analysis of public de-identified data | ☐ |
| B4 | Author Contributions (CRediT) — fill names | ☐ |
| B5 | Conflict of Interest — confirm by all authors | ☐ |
| B6 | Funding — fill or "no specific grant" statement | ☐ |
| B7 | AI Disclosure Statement — present | ☐ |
| B8 | References — Vancouver numbered, 16, all verified real (nature-ref-verifier, commit cf27614) | ☐ |
| B9 | Prognostic claims — "suggestive"/"borderline" wording + power/surrogate-endpoint caveats in Limitations | ☐ |

## C. Figures & tables

| # | Item | Spec | Done |
|---|------|------|------|
| C1 | Figures 1–5 + Supplementary Figure S3 exported | RGB; ≥300 dpi (color), ≥600 (mixed), ≥1000 (line); font ≥7 pt Arial | ☐ |
| C2 | Figure width | ≤183 mm (full column) | ☐ |
| C3 | Figure legends complete, self-contained | — | ☐ |
| C4 | Tables 1–3 + Supplementary Tables 1–3 | Numbered, cited in text | ☐ |
| C5 | Supplementary Methods (S1/S2: software versions, seeds, CMap screen) | — | ☐ |

## D. Files to upload (compile before submission)

- [ ] Manuscript (docx or LaTeX per journal preference) — with title page
- [ ] Cover letter (`docs/cover-letter.md`, fields filled)
- [ ] Figure files 1–5 + S3 (individual files, named per journal convention)
- [ ] Tables as separate files if required
- [ ] Supplementary files (Methods S1/S2, Tables S1–S3, Figure S3)
- [ ] Optional: suggested reviewers (3–5, no COI with authors)

## E. Reproducibility & archiving

| # | Item | Status | Done |
|---|------|--------|------|
| E1 | GitHub repository public? | **Yes — already public** (verified 2026-10-18); URL `github.com/PaperValidation/NF-PitNET-Methylation-Validation` is live for reviewers | ☑ |
| E2 | Zenodo DOI 10.5281/zenodo.21840024 | Restricted — **publish on acceptance**; creators placeholder "Tidus to be replaced" → replace with author list | ☐ |
| E3 | Raw data availability | Public (GEO/ArrayExpress/Mendeley) — no upload needed | ☐ |
| E4 | Analysis scripts committed | `scripts/` — path-canonicalized via `_paths.py`/`paths.R` (commit b5d7fd5 + follow-ups) | ☐ |

## F. Sign-off gates (require real people — cannot be automated)

- [ ] **AUTHOR_ACTION #1**: Full author list + affiliations + ORCIDs (manuscript line 18; Zenodo creators)
- [ ] **AUTHOR_ACTION #2**: Designate corresponding author (line 22; cover letter)
- [ ] **AUTHOR_ACTION #3**: CRediT role assignments (lines 360–364)
- [ ] **AUTHOR_ACTION #4**: COI confirmation by all authors (line 368)
- [ ] **AUTHOR_ACTION #5**: Funding statement (line 372)
- [ ] All authors have read and approved the final manuscript

> See `docs/author-action-checklist.md` for the exact per-line items.
