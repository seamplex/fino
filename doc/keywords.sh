UPPER1=`grep strcasecmp ../src/parser.c ../wasora/src/parser.c ../wasora/src/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[A-Z]/'`
UPPER2=`grep keywords   ../src/parser.c ../wasora/src/parser.c ../wasora/src/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[A-Z]/'`
LOWER=`grep strcasecmp  ../src/parser.c ../wasora/src/parser.c ../wasora/src/mesh/parser.c | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | mawk '/^[a-z]/' | sort`
VARS=`grep variable     ../src/init.c   ../wasora/src/init.c   ../wasora/src/mesh/init.c   | grep -v "computing" | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | sort`
FUNCS=`cat ../wasora/src/builtin.h                                                  | sed -r 's/[^"]*("[^"]*")?/ \1/g;s/" +"/\n/g;s/ *"//g'| mawk '$1 in p{next}{p[$1];print}' | sort`
UPPER=`echo $UPPER1 $UPPER2 | sort`
