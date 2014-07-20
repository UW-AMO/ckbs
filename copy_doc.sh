#! /bin/bash -eu
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# -------------------------------------------------------------------
# exit on any error
set -e
#
# check that this is the trunk
dir=`pwd | sed -e 's|.*/ckbs/||'`
if [ "$dir" != "trunk" ]
then
	echo "copy_doc.sh: can only be run in the trunk."
	echo "To change the documentation for a stable version,"
	echo "use ./build_doc.sh and commit changes in doc/*"
	exit 1
fi
#
echo "Determining current documentation version number"
version=`date +%F | sed -e 's|-||g' -e 's|^|0.|'`
echo "version = $version"
#
echo "Changing verison number in ckbs.omh"
sed -e "s/ckbs-[0-9.]\{10\}/ckbs-$version/" -i "" ckbs.omh 
#
echo "./build_doc.sh"
./build_doc.sh
#
echo "cd doc"
cd doc
#
if [ ! -e html ]
then
	echo "svn checkout \\"
	echo "	https://projects.coin-or.org/svn/CoinBazaar/html/ckbs html"
	svn checkout https://projects.coin-or.org/svn/CoinBazaar/html/ckbs html
else
	echo "svn update html"
	svn update html
fi
#
old_list=`ls html/* | sed -e 's|html/||'`
for file in $old_list
do
	if [ ! -e $file ]
	then
		echo "svn delete html/$file"
		svn delete html/$file
	fi
done
#
new_list=`ls * | sed -e 's|html:||'`
for file in $new_list
do
	if [ ! -e html/$file ]
	then
		echo "cp $file html/$file"
		cp $file html/$file
		#
		echo "svn add html/$file"
		svn add html/$file
	else
		echo "cp $file html/$file"
		cp $file html/$file
	fi
done
echo ""
echo "To see the differences use:"
echo "	svn diff doc/html"
echo "To commit the changes use:"
echo "	svn commit -m \"[html/ckbs] message\" doc/html"
echo "to commit the changes."
