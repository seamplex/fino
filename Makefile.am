plugin = fino

SUBDIRS = src
DIST_SUBDIRS = src wasora

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO
EXTRA_DIST = examples doc

TESTS = examples/test.sh

all-local:
	if [ -e src/${plugin} ]; then cp -u  src/${plugin} .; fi
	if [ -e ${plugin} ]; then cd examples; ln -sf ../${plugin}; cd;  fi
	if [ -e src/${plugin}-version.$(OBJEXT) ]; then rm -f src/${plugin}-version.$(OBJEXT);  fi

clean-local:
	rm -f src/${plugin}
	rm -f ${plugin} ${plugin}.so
