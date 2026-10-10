# 文献调研 skills（Claude Code 云端）

在这个仓库里开的 Claude Code 云端会话会自动加载下面两个 skill。

| skill | 用途 | 来源 |
|---|---|---|
| `paper-search-pro` | 同时检索 OpenAlex、Semantic Scholar、Crossref、arXiv、PubMed，AI 题目会标注 OpenReview 录用情况；子代理逐篇筛选，输出网页报告、BibTeX、RIS | [O0000-code/paper-search-pro](https://github.com/O0000-code/paper-search-pro) v2.5.0（commit `b551d74`），Apache-2.0 |
| `research-lit` | 找相关工作并整理成 related work，检索到的论文会逐篇核验，防止编造文献 | [wanshuiyin/Auto-claude-code-research-in-sleep](https://github.com/wanshuiyin/Auto-claude-code-research-in-sleep)（ARIS，commit `26b95cf`），MIT |

## 第一次使用前

### 1. 放开网络

云环境默认拦截学术网站。打开会话标题栏的云环境菜单 → Edit → Network access，选更宽的访问级别，或者在 Allowed domains 里加入下面这些域名（保留 Allow package managers 的勾选）：

```
api.openalex.org
api.semanticscholar.org
export.arxiv.org
arxiv.org
api.crossref.org
api2.openreview.net
openreview.net
doi.org
eutils.ncbi.nlm.nih.gov
```

### 2. 申请 key，填进环境变量

在同一个设置页添加环境变量，新开的会话才会生效。key 不要写进仓库，也不要贴到聊天里。

| 变量 | 是否必需 | 怎么拿 |
|---|---|---|
| `OPENALEX_API_KEY` | 必需（paper-search-pro 没有它会停下） | 在 https://openalex.org 注册，然后到 https://openalex.org/settings/api 生成，免费 |
| `RESEARCH_EMAIL` | 必需 | 填你的邮箱，不用申请。会作为联系邮箱发给 OpenAlex、Crossref、NCBI 的接口 |
| `SEMANTIC_SCHOLAR_API_KEY` | 推荐 | https://www.semanticscholar.org/product/api 申请，免费，审核可能要几天；没有也能用，只是限速更严 |

### 3. 连接 alphaXiv（推荐）

在 claude.ai 的「设置 → 连接器」里搜 alphaXiv 并登录，就能在 arXiv 全文里检索、对 PDF 提问。

## 用法

- 直接说需求，会触发 paper-search-pro：「帮我找 2023 年以来用 neural operator 求解偏微分方程的论文」。默认是标准深度（约 10–17 分钟，60–180 篇）；说「快速扫一下」或「做综述，要全面」可以换深度。结果在 `paper-search-results/<检索编号>/report.html`。
- 整理相关工作用 `/research-lit`：
  `/research-lit "equivariant graph neural networks for molecular dynamics" — sources: all, semantic-scholar, openalex`
  不加 `— sources:` 时用网页搜索、arXiv 接口和本地 PDF（放在 `papers/` 或 `literature/` 目录）；Semantic Scholar 和 OpenAlex 要在 `— sources:` 里写明才会查。

## 文件说明

- `.claude/skills/paper-search-pro/`、`.claude/skills/research-lit/`：两个 skill 本体。
- `.claude/skills/shared-references/`：research-lit 引用的规范文档（来自 ARIS）。
- `.aris/tools/`：research-lit 的辅助脚本，负责 arXiv、Semantic Scholar、OpenAlex 查询和引用核验（来自 ARIS）。
- `.claude/hooks/session-start.sh`：每次开云端会话时安装 Python 依赖，并用上面的环境变量生成 `~/.paper-search-pro/config.yaml`；还会报告哪些 key 缺失、哪些网站被网络策略拦住。

## 更新

两个 skill 都是从上游复制过来的。要更新时跟 Claude 说「把 paper-search-pro（或 research-lit）更新到上游最新版」即可。
