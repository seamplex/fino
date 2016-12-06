if [ -e ./runtest.sh ]; then
 . runtest.sh
elif [ -e examples/runtest.sh ]; then
 cd examples
 . runtest.sh
else
 echo wrong PWD
 exit 77
fi
