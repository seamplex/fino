plugin = fino

SUBDIRS = src
DIST_SUBDIRS = src wasora

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO
EXTRA_DIST = examples doc

TESTS = examples/test.sh

all-local:
	if [ -e src/.libs/lib${plugin}.so.0.0.0 ]; then cp -u  src/.libs/lib${plugin}.so.0.0.0 ./${plugin}.so; fi
	if [ -e src/wasora/${plugin}$(EXEEXT) ]; then cp -u  src/wasora/${plugin}$(EXEEXT) .; fi
	if [ -e ${plugin}$(EXEEXT) ]; then cd examples; ln -sf ../${plugin}$(EXEEXT); cd;  fi

clean-local:
	rm -f src/wasora/${plugin}$(EXEEXT)
	rm -f ${plugin} ${plugin}.so
