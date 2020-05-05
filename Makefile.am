plugin = fino

SUBDIRS = src doc
DIST_SUBDIRS = src wasora doc

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO \
                doc/fino.texi doc/fino.info doc/fino.svg doc/fino.xml doc/reference.txt

dist_man_MANS = doc/fino.1

EXTRA_DIST = examples locateruntest.sh \
             src/version.sh

TESTS = examples/test-tensile.sh \
        examples/test-piston.sh \
        examples/test-wire.sh \
        examples/test-heated-cylinder.sh \
        examples/test-cantilever.sh
        
all-local:
	if [ -e src/${plugin} ]; then cp  src/${plugin} .; fi
	if [ -e ${plugin} ]; then cd examples; ln -sf ../${plugin}; cd;  fi
# 	if [ -e src/${plugin}-version.$(OBJEXT) ]; then rm -f src/${plugin}-version.$(OBJEXT);  fi

clean-local:
	rm -f src/${plugin}
	rm -f ${plugin} ${plugin}.so
