\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{abs}[2021]
% article inheritance
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
% packages
\RequirePackage[margin=25mm,left=45mm]{geometry}
\RequirePackage{authblk}
% spacing
\setlength{\parskip}{1ex}
\setlength{\parindent}{0ex}
\def\@parboxrestore{\setlength{\parskip}{1ex}}
\frenchspacing
\emergencystretch=2em
% meta-data
\renewcommand{\title}[2][]{\gdef\@maxtitle{#1}\gdef\@title{#2}}
\newcommand{\conference}[2][]{\gdef\@conference{#2 #1}}
\renewcommand{\abstract}[2][]{\gdef\@maxwords{#1}\gdef\@abstract{#2}}
\newcommand{\keywords}[1]{\gdef\@keywords{\ignorespaces#1}}
\newcommand{\float}[3][Figure]{\gdef\@floatlab{#1}\gdef\@floatcap{#2}\gdef\@float{#3}}
% helpers
\renewcommand{\section}[1]{\par\MakeUppercase{#1.}}
\newcommand{\tstyle}[1]{\textbf{#1}}
\newcommand{\tpart}[3][:]{%
  \makebox[0pt][r]{\tstyle{#2#1}\hspace{1mm}}#3\par\smallskip}
\def\sep{\unskip\hskip.4em\raisebox{.5ex}{\textbf.}\hskip.4em}
% maketitle
\renewcommand{\maketitle}{%
  \immediate\write18{detex \@abstract.tex | wc -w > .tmp/abswords.tex}
  \immediate\write18{detex \@title.tex | wc -c > .tmp/titchars.tex}
  \thispagestyle{empty}
  \tpart{Title}{\input\@title\ignorespaces}
  \tpart{Authors}{\@author}
  \tpart{Conference}{\@conference}
  \tpart{Submitted}{\@date}
  \tpart{Abstract}{\vspace{-\parskip}\vspace{-\baselineskip}%
    \input\@abstract\ignorespaces\par\bigskip}
  \tpart{Limits}{%
    \input{.tmp/titchars}\unskip\,/\,\@maxtitle\ \ title characters\par
    \input{.tmp/abswords}\unskip\,/\,\@maxwords\ \ abstract words}
  \tpart{Keywords}{\@keywords}
  \pagebreak
  \tpart{\@floatlab}{\@floatcap\par
    \hskip-10mm\parbox{\linewidth}{\@float}}}
