# generates READMEs and INSTALLs if different formats from the markdown source

if [ -z "`which pandoc`" ]; then
 echo "pandoc is not installed"
 exit 1
fi

pandoc README.md -t plain -o README
pandoc README.md -o doc/README.pdf -V geometry:margin=3cm # --default-image-extension=pdf
pandoc README.md -o doc/README.html -s                    # --default-image-extension=svg
pandoc INSTALL.md -t plain -o INSTALL
pandoc INSTALL.md -o doc/INSTALL.pdf -V geometry:margin=3cm
pandoc INSTALL.md -o doc/INSTALL.html -s
pandoc ${WASORA_PATH}/PLUGINS.md -t plain -o PLUGINS
pandoc ${WASORA_PATH}/PLUGINS.md -o doc/PLUGINS.pdf
pandoc ${WASORA_PATH}/PLUGINS.md -o doc/PLUGINS.html -s -S

# cd doc
# pandoc README.md -o README.pdf -V geometry:margin=3cm # --default-image-extension=pdf
# pandoc README.md -o README.html -s                    # --default-image-extension=svg
# cd ..
