#!/bin/sh
# 
# $ ./syntax-tex.sh > ../syntax.tex
# 

. ./keywords.sh

cat << EOF
\definecolor{was_fondo}{rgb}{0.95,0.95,0.90}
\definecolor{was_keyword1}{rgb}{0.0,0.0,0.0}
\definecolor{was_keyword2}{rgb}{0.0,0.2,0.0}
\definecolor{was_variable}{rgb}{0.5,0.2,0.2}
\definecolor{was_function}{rgb}{0.2,0.5,0.2}
\definecolor{was_comment}{rgb}{0.5,0.5,0.5}

\definecolor{bash_fondo}{rgb}{0.97,0.97,0.97}
\definecolor{bash_keyword1}{rgb}{0.9,0.9,0.9}
\definecolor{bash_keyword2}{rgb}{0.7,0.7,0.7}
\definecolor{bash_comment}{rgb}{0.5,0.5,0.5}
EOF

# primary keywords
echo "\lstdefinelanguage{fino}{"
echo "morekeywords={"
for kw in $UPPER; do
  echo "      $kw,"
done
echo "},"

# secondary keywords (TO-DO)
echo "morekeywords={[2]"
echo "},"

# special variables
echo "morekeywords={[3]"
for kw in $VARS; do
  echo "      $kw,"
  echo -n "      $kw"
  echo "_0,"
done
echo "},"

# functions
echo "morekeywords={[4]"
for kw in $FUNCS; do
  echo "      $kw,"
done
echo "},"

cat << EOF
sensitive=true,
morecomment=[l]{\#},
morestring=[b]\",
}

\newcommand{\MyHookSign}{\hbox{\ensuremath{\hookleftarrow}}}

\lstset{
  language=fino,
  basicstyle=\footnotesize,
  commentstyle={\color{was_comment}\textit},
  keywordstyle=[1]{\color{was_keyword1}\ttfamily\textbf},
  keywordstyle=[2]{\color{was_keyword2}\ttfamily\textbf},
  keywordstyle=[3]{\color{was_variable}\textit},
  keywordstyle=[4]{\color{was_function}\textbf},
  backgroundcolor=\color{was_fondo},
  breaklines=true,
  prebreak={\space\MyHookSign},
  xleftmargin=0.8cm,
  xrightmargin=0.5cm,
  frame=single,
  framesep=0.5cm
}


\lstdefinestyle{fino}{
  language=fino,
  basicstyle=\footnotesize,
  commentstyle={\color{was_comment}\textit},
  keywordstyle=[1]{\color{was_keyword1}\ttfamily\textbf},
  keywordstyle=[2]{\color{was_keyword2}\ttfamily\textbf},
  keywordstyle=[3]{\color{was_variable}\textit},
  keywordstyle=[4]{\color{was_function}\textbf},
  backgroundcolor=\color{was_fondo},
  breaklines=true,
  prebreak={\space\MyHookSign},
  xleftmargin=0.5cm,
  xrightmargin=0.5cm,
  framesep=0.5cm,
  frame=single,
}

\lstdefinestyle{bash}{
  language=,
  basicstyle=\footnotesize\ttfamily,
  backgroundcolor=\color{bash_fondo},
  breaklines=true,
  prebreak={\space\MyHookSign},
  xleftmargin=0.2cm,
  xrightmargin=0.2cm,
  framesep=0.2cm,
  frame=single,
}
EOF
