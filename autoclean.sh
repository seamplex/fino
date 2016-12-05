#!/bin/sh
# 
# Execute this script to clean the directory and leave it as
# a fresh mercurial branch
#
m4 m4/bootstrap.m4 - << EOF | sh
plugin=`grep plugin= autogen.sh | head -n1 | cut -c8-`
PLUGIN_AUTOCLEAN
EOF
