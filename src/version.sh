rm -f version.h
if test -e ../.hg -a ! -z "`which hg`"; then
  majorversion=`hg log -r tip --template='{latesttag}'`
  minorversion=`hg log -r tip --template="{latesttagdistance}"`
  hg log -r tip --template="\
#define PLUGIN_VCS_BRANCH    \"{branch}\"\n\
#define PLUGIN_VCS_MAJOR     \"${majorversion}\"\n\
#define PLUGIN_VCS_MINOR     \"${minorversion}\"\n\
#define PLUGIN_VCS_REVID     \"{node}\"\n\
#define PLUGIN_VCS_SHORTID   \"{node|short}\"\n\
#define PLUGIN_VCS_REVNO     {rev}\n\
#define PLUGIN_VCS_DATE      \"{date|isodate}\"\n" > version-vcs.h
fi

cat version-vcs.h >> version.h

if [ -e version-conf.h ]; then
  cat version-conf.h >> version.h
fi

cat << EOF >> version.h
#define COMPILATION_DATE     "`date +'%Y-%m-%d %H:%M:%S'`"
#define COMPILATION_USERNAME "`whoami | sed s/\\\\\\\\//`"
#define COMPILATION_HOSTNAME "`hostname`"
#define PLUGIN_DATE          "`stat -c %y *.c | sort -r | head -n1 | cut -b-19`"
#define PLUGIN_HEADERMD5     "`md5sum wasora/wasora.h | cut -c-32`"
EOF
