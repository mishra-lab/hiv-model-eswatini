\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{poster}
[2021/09/01 v0.0.0 custom poster class based on tikzposter]
% spacing & tikzposter stuff
\newlength{\gap}\setlength{\gap}{1cm}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{tikzposter}}
\PassOptionsToClass{
  a0paper,
  margin=0cm,
  innermargin=\gap,
  colspace=\gap,
  subcolspace=\gap,
  blockverticalspace=\gap,
}{tikzposter}
\ProcessOptions\relax
\LoadClass{tikzposter}
\geometry{paperwidth=1200mm,paperheight=600mm}
\setlength{\TP@visibletextwidth}{\textwidth-2\TP@innermargin}
\setlength{\TP@visibletextheight}{\textheight-2\TP@innermargin}
\tikzposterlatexaffectionproofoff
\useblockstyle[
  linewidth=.2\gap,
  titleinnersep=.5\gap,
  bodyinnersep=\gap,
  roundedcorners=.5\gap,
]{Default}
% maketitle + logos
\let\omaketitle\maketitle
\AtBeginDocument{\renewcommand{\maketitle}{%
  \omaketitle[%
    innersep=0mm,
    width=\dimexpr\paperwidth-2\gap,
    linewidth=.2\gap,
    roundedcorners=.5\gap,
    titletotopverticalspace=\gap,
    titletoblockverticalspace=\gap,
    titlegraphictotitledistance=0cm]
  \begin{pgfonlayer}{notelayer}
    \node[anchor=north west,align=left] at
      (\dimexpr+1.5\gap-.5\paperwidth,\dimexpr-1.5\gap+.5\paperheight){\logo@left};
    \node[anchor=north east,align=right] at
      (\dimexpr-1.5\gap+.5\paperwidth,\dimexpr-1.5\gap+.5\paperheight){\logo@right};
  \end{pgfonlayer}}}
\newcommand{\logoleft}[1]{\gdef\logo@left{#1}}\logoleft{}
\newcommand{\logoright}[1]{\gdef\logo@right{#1}}\logoright{}
\newcommand{\titlestyle}[1]{\gdef\title@style{#1}}\titlestyle{\Huge\scshape}
\settitle{\centering\vbox{
  \color{titlefgcolor}
  \vspace*{2ex}
  {\title@style\selectfont\@title\par}
  \vspace*{1ex}
  {\LARGE\@author\par}
  \vspace*{1ex}
  {\large\@institute}
  \vspace*{2ex}
}}
\newcommand{\af}[1]{\textsuperscript{\,#1}}
\newcommand{\affil}[2]{\textsuperscript{#1\,}#2\hskip1ex}
% colors
\definecolor{C1}{HTML}{002244}
\definecolor{C2}{HTML}{CC0033}
\definecolor{C3}{HTML}{0099CC}
\definecolor{C4}{HTML}{ff9900}
\definecolor{C5}{HTML}{990099}
\usecolorstyle[
  colorOne=C1!10!white,
  colorTwo=C1,
  colorThree=C2,
]{Default}
\colorlet{titlebgcolor}{C1}
\colorlet{titlefgcolor}{white}
\RequirePackage[
  colorlinks=true,
  urlcolor=C3,
  linkcolor=C3,
]{hyperref}
% bib
\RequirePackage[
  backend=biber,
  style=numeric,
  sorting=none,
  natbib=true,
  url=false,
  doi=true,
  isbn=false,
  locallabelwidth,
]{biblatex}
\renewcommand{\bibfont}{\footnotesize}
\renewcommand\@biblabel[1]{#1}
\let\oprintbibliography\printbibliography
\def\printbibliography{%
  \oprintbibliography[heading=none]}
% figures
\newcounter{fig}
\newenvironment{fig}[2][]{
  \begin{minipage}[#1]{#2}\small}{\end{minipage}}
\renewcommand{\caption}[1]{%
  \refstepcounter{fig}%
  {\small\textbf{Figure~\thefig.} #1\par}}
\def\HyperFirstAtBeginDocument#1{#1} % BUG
