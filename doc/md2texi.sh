if [ -z "`which pandoc`" ]; then 
 echo "error: pandoc not installed"
 exit 1
fi

pandoc help.md -t plain > help.txt
m4 header.m4 reference-manual.m4 > reference-manual.md
m4 header.m4 fino.1.md | pandoc -s -t man -o fino.1
m4 header.m4 fino.md | pandoc --toc --template template.texi -o fino.texi

