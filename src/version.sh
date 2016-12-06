rm -f version.h

if test -e ../.git -a ! -z "`which git`"; then
 version=`git describe | sed 's/-/./'`
 echo "[[define]](wasoraversion, ${version})[[dnl]]" > version.m4

 branch=$(git symbolic-ref HEAD | sed -e 's,.*/\(.*\),\1,')
 date=`git log --pretty=format:"%ad" | head -n1`
 cat << EOF > version-vcs.h
#define PLUGIN_VCS_BRANCH    "${branch}"
#define PLUGIN_VCS_VERSION   "${version}"
#define PLUGIN_VCS_DATE      "${date}"
#define PLUGIN_VCS_CLEAN     `git status --porcelain | wc -l`
EOF
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
#define PLUGIN_HEADERMD5     "`md5sum ../wasora/src/wasora.h | cut -c-32`"
EOF
