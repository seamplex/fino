if [ -z "`which pandoc`" ]; then 
 echo "error: pandoc not installed"
 exit 1
fi

# m4 reference.m4 > reference.md
pandoc help.md -t plain > help.txt

m4 header.m4 fino.md | pandoc --toc --template template.texi -o fino.texi

m4 header.m4 fino.md > tmp
# pandoc tmp.md --toc --template template.texi -o blackjack.texi
