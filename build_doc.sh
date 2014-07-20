#! /bin/bash -e
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# -------------------------------------------------------------------
#
count=`ls | grep omhelp- | wc -w | sed -e 's| ||g'`
if [ "$count" = 1 ]
then
	dir=`ls | grep omhelp-`
else
	dir="do_not_use_current_omhelp_directory"
fi
if [ -e $dir/src/omhelp ]
then
	echo "using existing $dir/src/omhelp"
else
	if [ ! -e OMhelp.unix.tar.gz ]
	then
		check=`which wget | sed -e 's|.*/||'`
		if [ "$check" != "wget" ]
		then
			echo "build_doc.sh: cannot find wget command"
			exit 1
		fi
		web_page="http://www.seanet.com/~bradbell"
		echo "wget $web_page/OMhelp.unix.tar.gz"
		wget "$web_page/OMhelp.unix.tar.gz"
	fi
	if [ -e "omhelp-*" ]
	then
		echo "rm -rf omhelp-*"
		rm -rf omhelp-*
	fi
	echo "tar -xvzf OMhelp.unix.tar.gz"
	tar -xvzf OMhelp.unix.tar.gz
	#
	echo "cd omhelp-*"
	cd omhelp-*
	#
	echo "./configure --prefix=$HOME"
	./configure --prefix=$HOME
	#
	echo "make"
	make
	#
	echo "cd .."
	cd ..
fi
count=`ls | grep omhelp- | wc -w | sed -e 's| ||g'`
if [ "$count" != 1 ]
then
	echo "build_doc.sh: unknown error, giving up"
	exit 1
fi
dir=`ls | grep omhelp-`
if [ ! -e doc ]
then
	echo "mkdir doc"
	mkdir doc
fi
#
echo "cd doc"
cd doc
#
# build the documentation
home_page="https://projects.coin-or.org/CoinBazaar/wiki/Projects/ckbs"
cmd="../$dir/src/omhelp ../ckbs.omh -omhelp_dir ../$dir/OMhelp -image_link $home_page"
log="../omhelp.log"
for arg1 in "" -xml
do
	for arg2 in -noframe -printable
	do
		echo "omhelp ckbs.omh $arg1 $arg2 >& omhelp.log"
		if ! $cmd $arg1 $arg2 -debug >& ../omhelp.log
		then
			echo "build_doc.sh: omhelp error, see omhelp.log"
			exit 1
		fi
		if grep "^OMhelp Warning:" ../omhelp.log
		then
			echo "build_doc.sh: omhelp warning, see omhelp.log"
			exit 1
		fi
	done
done
