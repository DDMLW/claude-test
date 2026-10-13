"""
build.py —— 由 manuscript.md 生成 Chaos, Solitons & Fractals（Elsevier elsarticle）投稿格式的 LaTeX 并编译
================================================================================================
manuscript.md 是正文的唯一源文件；本脚本只做格式转换，不改内容：
  1. 按一级标题拆分：Abstract / Highlights / Keywords / 正文各节 / Appendix A、B / Figure captions /
     Supplementary figure captions / 后置节（CRediT、Declaration、Data availability、Acknowledgments）；
  2. 每段用 pandoc（--natbib）转成 LaTeX：[@a; @b] → \\citep{a,b}；删除 HTML 注释与 [C#] 出处标注；去掉标题前的编号；
  3. 图注取自 “**Fig. N. 标题.** 正文” 段落，按 FIGS 放在对应节之后（浮动体）；SI 图放在附录之后；
  4. 写 latex/main.tex（elsarticle preprint 12pt、数字引用 elsarticle-num）、复制 refs.bib，pdflatex + bibtex 编译 → latex/main.pdf。
需要：pandoc；编译还需要 pdflatex、bibtex 与 elsarticle（TeX Live 的 texlive-publishers）。
运行：python build.py            # 生成并编译
      python build.py --no_pdf   # 只生成 latex/main.tex 与 latex/body_*.tex
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "latex")
# 正文第几节之后放哪些图: {节号: [(图号, [文件名, ...]), ...]}
FIGS = {"3": [("1", ["symmetry_transfer.png"])],
        "4": [("2", ["atlas_fields.png"]), ("3", ["scaling_local_slopes.png"])],
        "5": [("4", ["window_mechanism.png"])],
        "6": [("5", ["dynamic_window.png"])],
        "7": [("6", ["twin_vortex_fields.png", "twin_vortex.png"])]}
SI_FIGS = [("S1", ["e3c_radial_vs_2d.png"]), ("S2", ["twin_vortex_seed.png"])]
GRAPHICS = ["../../results/analysis/hq/", "../../results/E3_stability/hq/"]
BACK = ("CRediT authorship contribution statement", "Declaration of competing interest", "Data availability",
        "Acknowledgments")

PREAMBLE = r"""\documentclass[preprint,12pt]{elsarticle}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{longtable,booktabs,array,calc}
\usepackage{lineno}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\graphicspath{%(graphics)s}
\journal{Chaos, Solitons \& Fractals}
\begin{document}
\begin{frontmatter}
\title{%(title)s}
\author{%(author)s}
\begin{abstract}
%(abstract)s
\end{abstract}
\begin{highlights}
%(highlights)s
\end{highlights}
\begin{keyword}
%(keywords)s
\end{keyword}
\end{frontmatter}
\linenumbers
"""


def pandoc(md, inline=False):
    """Markdown 片段 → LaTeX（natbib 引用）。"""
    r = subprocess.run(["pandoc", "-f", "markdown", "-t", "latex", "--natbib", "--wrap=preserve"],
                       input=md, capture_output=True, text=True, check=True)
    tex = r.stdout.strip()
    if inline:          # 只有一个段落时去掉首尾空行
        tex = tex.replace("\n\n", " ")
    return tex


def clean(md):
    md = re.sub(r"<!--.*?-->", "", md, flags=re.S)
    md = re.sub(r" ?\[C\d+(?:, ?C\d+)*\]", "", md)          # 出处标注 [C8]
    return md.strip() + "\n"


def split_sections(text):
    """[(一级标题, 内容)]，内容不含标题行。"""
    parts = re.split(r"^# (.+)$", text, flags=re.M)
    return [(parts[i].strip(), parts[i + 1]) for i in range(1, len(parts), 2)]


def captions(md):
    """“**Fig. N. 标题.** 正文” → {N: LaTeX 图注（不含 “Fig. N.”）}。"""
    out = {}
    for blk in re.split(r"\n(?=\*\*Fig\. )", "\n" + md):
        m = re.match(r"\s*\*\*Fig\. (S?\d+)\. (.*?)\*\*(.*)", blk, flags=re.S)
        if m:
            out[m.group(1)] = pandoc(clean(f"**{m.group(2).strip()}**{m.group(3).rstrip()}"), inline=True)
    return out


def find_graphic(f):
    """GRAPHICS 中第一个存在该文件的目录 (相对 latex/); 没有时 None。"""
    for g in GRAPHICS:
        if os.path.exists(os.path.normpath(os.path.join(OUT, g, f))):
            return g
    return None


def figure(num, files, cap):
    """图的浮动体; 缺少的图片用占位框代替 (例如 E3c 尚未在 hq 档运行时), 不让编译失败。"""
    imgs = "\n\\par\\medskip\n".join(
        (f"\\includegraphics[width=0.98\\textwidth,height=0.72\\textheight,keepaspectratio]{{{f}}}" if len(files) == 1 else
         f"\\includegraphics[width=0.98\\textwidth,height=0.4\\textheight,keepaspectratio]{{{f}}}") if find_graphic(f) else
        f"\\fbox{{\\parbox{{0.9\\textwidth}}{{\\centering missing figure: \\texttt{{{f.replace('_', chr(92) + '_')}}}}}}}"
        for f in files)
    label = f"fig{num}"
    return (f"\\begin{{figure}}[tbp]\n\\centering\n{imgs}\n\\caption{{{cap}}}\n\\label{{{label}}}\n\\end{{figure}}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--no_pdf", action="store_true", help="只生成 tex，不编译")
    a = p.parse_args()
    src = open(os.path.join(HERE, "manuscript.md"), encoding="utf-8").read()
    head = re.match(r"---\n(.*?)\n---\n", src, flags=re.S)
    title = re.search(r'^title:\s*"(.*)"', head.group(1), flags=re.M).group(1)
    author = re.search(r'^author:\s*"(.*)"', head.group(1), flags=re.M).group(1)
    secs = split_sections(src[head.end():])
    get = {t: c for t, c in secs}
    abstract = pandoc(clean(get["Abstract"]))
    highlights = "\n".join("\\item " + pandoc(clean(ln[2:]), inline=True)
                           for ln in get["Highlights"].strip().splitlines() if ln.startswith("* "))
    keywords = " \\sep ".join(k.strip() for k in clean(get["Keywords"]).strip().split(";"))
    caps = {**captions(get.get("Figure captions", "")), **captions(get.get("Supplementary figure captions", ""))}
    os.makedirs(OUT, exist_ok=True)
    body, appendix, back = [], [], []
    for t, c in secs:
        m = re.match(r"(\d+)\. (.+)", t)
        if m:                                            # 正文
            md = clean(f"# {m.group(2)}\n\n" + re.sub(r"^## \d+\.\d+ ", "## ", c, flags=re.M))
            body.append(pandoc(md))
            for num, files in FIGS.get(m.group(1), []):
                body.append(figure(num, files, caps.get(num, "")))
        elif t.startswith("Appendix"):
            name = re.sub(r"^Appendix [A-Z]\.\s*", "", t)
            appendix.append(pandoc(clean(f"# {name}\n\n{c}")))
        elif t in BACK:
            back.append(pandoc(clean(f"# {t}\n\n{c}")).replace("\\section{", "\\section*{"))
    si = ("\\section{Supplementary figures}\n\\renewcommand{\\thefigure}{S\\arabic{figure}}\\setcounter{figure}{0}\n"
          + "".join(figure(n, f, caps.get(n, "")) for n, f in SI_FIGS))
    parts = {"body": "\n\n".join(body), "back": "\n\n".join(back),
             "appendix": "\\appendix\n" + "\n\n".join(appendix) + "\n\n" + si}
    for k, v in parts.items():
        with open(os.path.join(OUT, f"body_{k}.tex"), "w", encoding="utf-8") as fh:
            fh.write(v + "\n")
    main_tex = PREAMBLE % dict(graphics="".join("{%s}" % g for g in GRAPHICS), title=title, author=author,
                               abstract=abstract, highlights=highlights, keywords=keywords)
    main_tex += ("\\input{body_body}\n\\input{body_back}\n\\input{body_appendix}\n"
                 "\\bibliographystyle{elsarticle-num}\n\\bibliography{refs}\n\\end{document}\n")
    with open(os.path.join(OUT, "main.tex"), "w", encoding="utf-8") as fh:
        fh.write(main_tex)
    shutil.copy(os.path.join(HERE, "refs.bib"), os.path.join(OUT, "refs.bib"))
    print(f"written {OUT}/main.tex (+ body_body/back/appendix.tex, refs.bib)")
    if a.no_pdf:
        return
    log = []
    for cmd in (["pdflatex", "-interaction=nonstopmode", "main"], ["bibtex", "main"],
                ["pdflatex", "-interaction=nonstopmode", "main"], ["pdflatex", "-interaction=nonstopmode", "main"]):
        r = subprocess.run(cmd, cwd=OUT, capture_output=True, text=True)
        log.append(r.stdout[-3000:])
        if r.returncode != 0 and cmd[0] == "pdflatex":
            print(r.stdout[-4000:])
            sys.exit(f"{' '.join(cmd)} failed")
    txt = open(os.path.join(OUT, "main.log"), encoding="latin-1").read()
    warn = sorted(set(re.findall(r"(Citation `[^']+' .*?undefined|Reference `[^']+' .*?undefined)", txt)))
    print(f"latex/main.pdf built; {len(warn)} undefined citations/references")
    for w in warn[:20]:
        print("  ", w)


if __name__ == "__main__":
    main()
