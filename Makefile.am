plugin = fino

SUBDIRS = src
DIST_SUBDIRS = src wasora

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO
EXTRA_DIST = examples doc

TESTS = examples/test.sh

all-local:
# 	if [ -e src/.libs/lib${plugin}.so.0.0.0 ]; then cp -u  src/.libs/lib${plugin}.so.0.0.0 ./${plugin}.so; fi
	if [ -e src/${plugin} ]; then cp -u  src/${plugin} .; fi
	if [ -e ${plugin} ]; then cd examples; ln -sf ../${plugin}; cd;  fi

clean-local:
	rm -f src/${plugin}
	rm -f ${plugin} ${plugin}.so
