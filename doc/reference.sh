touch reference-toc.md
m4 reference.m4 > reference.md
pandoc reference.md --toc --template=toc.template -o reference-toc.md
m4 reference.m4 > reference.md
