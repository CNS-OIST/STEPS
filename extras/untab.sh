#!/bin/sh

# $Id$

expand -t 4 $1 > $1.tmp
mv $1.tmp $1

# END
