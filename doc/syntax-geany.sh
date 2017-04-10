#!/bin/sh
# take the output of this script as a geany syntax definition file for fino
#
# ./syntax-geany.sh > $HOME/.config/geany/filedefs/filetypes.Fino.conf
#
# remember to add the line
# Fino=*.fin;
# to filetype_extensions.conf (tools -> configuration files -> filetype_extensions.conf)


. ./keywords.sh

cat << EOF
# For complete documentation of this file, please see Geany's main documentation
[styling]
# Edit these in the colorscheme .conf file instead
default=default
comment=comment
kword=keyword_1
operator=operator
basekword=keyword_2
otherkword=keyword_3
number=number_1
string=string_1
string2=string_2
identifier=identifier
infix=function
infixeol=function
package=function
package_other=keyword_2

[keywords]
# all items must be in one line
primary=${UPPER}
package_other=`echo ${VARS} | xargs`
package=`echo ${FUNCS} | xargs`

[settings]
# default extension used when saving files
extension=was
lexer_filetype=R

# the following characters are these which a "word" can contains, see documentation
#wordchars=_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789

# single comments, like # in this file
comment_single=#
# multiline comments
#comment_open=
#comment_close=

# set to false if a comment character/string should start at column 0 of a line, true uses any
# indentation of the line, e.g. setting to true causes the following on pressing CTRL+d
   #command_example();
# setting to false would generate this
#  command_example();
# This setting works only for single line comments
comment_use_indent=true

# context action command (please see Geany's main documentation for details)
context_action_cmd=

[indentation]
#width=4
# 0 is spaces, 1 is tabs, 2 is tab & spaces
type=0

[build_settings]
run_cmd=fino %f
EOF

# fortran
# cat << EOF
# # For complete documentation of this file, please see Geany's main documentation
# [styling]
# # Edit these in the colorscheme .conf file instead
# default=default
# comment=comment
# number=number_1
# string=string_1
# operator=operator
# identifier=identifier_1
# string2=string_2
# word=keyword_1
# word2=keyword_2
# word3=keyword_3
# preprocessor=preprocessor
# operator2=operator
# continuation=default
# stringeol=string_eol
# label=type
# 
# [keywords]
# # all items must be in one line
# primary=
# intrinsic_functions=
# user_functions=
# 
# [settings]
# # default extension used when saving files
# extension=was
# lexer_filetype=F77
# 
# # the following characters are these which a "word" can contains, see documentation
# #wordchars=_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
# 
# # single comments, like # in this file
# comment_single=#
# # multiline comments
# #comment_open=
# #comment_close=
# 
# # set to false if a comment character/string should start at column 0 of a line, true uses any
# # indentation of the line, e.g. setting to true causes the following on pressing CTRL+d
#         #command_example();
# # setting to false would generate this
# #       command_example();
# # This setting works only for single line comments
# comment_use_indent=false
# 
# # context action command (please see Geany's main documentation for details)
# context_action_cmd=
# 
# [indentation]
# #width=4
# # 0 is spaces, 1 is tabs, 2 is tab & spaces
# type=0
# 
# [build_settings]
# run_cmd=fion %f
# EOF
# 

# python
# cat << EOF
# # For complete documentation of this file, please see Geany's main documentation
# [styling]
# # Edit these in the colorscheme .conf file instead
# default=default
# commentline=comment_line
# number=number_1
# string=string_1
# character=character
# word=keyword_1
# triple=string_2
# tripledouble=string_2
# classname=type
# defname=function
# operator=operator
# identifier=function
# commentblock=comment
# stringeol=string_eol
# word2=keyword_2
# decorator=decorator
# 
# [keywords]
# # all items must be in one line
# primary=${UPPER}
# # identifier=`echo ${VARS} | xargs`
# # defname=`echo ${FUNCS} | xargs`
# 
# [lexer_properties]
# #fold.comment.python=1
# #fold.quotes.python=1
# 
# [settings]
# # default extension used when saving files
# extension=.was
# lexer_filetype=Python
# 
# 
# # the following characters are these which a "word" can contains, see documentation
# #wordchars=_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
# 
# # single comment char, like # in this file
# comment_single=#
# # multiline comments
# #comment_open="""
# #comment_close="""
# 
# # set to false if a comment character/string should start at column 0 of a line, true uses any
# # indentation of the line, e.g. setting to true causes the following on pressing CTRL+d
#         #command_example();
# # setting to false would generate this
# #       command_example();
# # This setting works only for single line comments
# comment_use_indent=true
# 
# # context action command (please see Geany's main documentation for details)
# context_action_cmd=
# 
# [indentation]
# #width=4
# # 0 is spaces, 1 is tabs, 2 is tab & spaces
# #type=0
# 
# [build_settings]
# # %f will be replaced by the complete filename
# # %e will be replaced by the filename without extension
# # (use only one of it at one time)
# compiler=python -m py_compile "%f"
# run_cmd=fino %f
# EOF


# # Sh
# cat << EOF
# # For complete documentation of this file, please see Geany's main documentation
# [styling]
# # Edit these in the colorscheme .conf file instead
# default=default
# commentline=comment_line
# number=number_1
# word=keyword_1
# string=string_1
# character=character
# operator=operator
# identifier=identifier_1
# backticks=backticks
# param=parameter
# scalar=scalar
# error=error
# here_delim=here_doc
# here_q=here_doc
# 
# [keywords]
# primary=${UPPER}
# identifier=`echo ${VARS} | xargs` `echo ${FUNCS} | xargs`
# 
# [settings]
# lexer_filetype=Sh
# # default extension used when saving files
# extension=.was
# 
# # the following characters are these which a "word" can contains, see documentation
# #wordchars=_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
# 
# # single comments, like # in this file
# comment_single=#
# # multiline comments
# #comment_open=
# #comment_close=
# 
# # set to false if a comment character/string should start a column 0 of a line, true uses any
# # indentation of the line, e.g. setting to true causes the following on pressing CTRL+d
#         #command_example();
# # setting to false would generate this
# #       command_example();
# # This setting works only for single line comments
# comment_use_indent=true
# 
# # context action command (please see Geany's main documentation for details)
# context_action_cmd=
# 
# [indentation]
# #width=4
# # 0 is spaces, 1 is tabs, 2 is tab & spaces
# #type=1
# 
# [build_settings]
# # %f will be replaced by the complete filename
# # %e will be replaced by the filename without extension
# # (use only one of it at one time)
# run_cmd="./%f"
# EOF
