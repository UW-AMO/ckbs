#! /bin/bash
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# ----------------------------------------------------------------------------
# exit on any error
set -e
#
# check initial working directory
dir=`pwd | sed -e 's|.*/||'`
if [ "$dir" != "trunk" ]
then
	echo "new_stable.sh: must execute this script in the trunk"
	exit 1
fi
#
echo "Getting current repository revision number"
rev_trunk=`svn info --revision HEAD | \
	grep '^Revision:' | sed -e 's|^Revision: *||'`
echo "rev_trunk = $rev_trunk"
#
echo "Getting current stable version number"
stable_version=`date +%F | sed -e 's/-//g' -e 's|^|0.|'`
echo "stable_version = $stable_version"
#
echo "Getting current documentation verison number"
doc_version=`grep '$section' ckbs.omh | \
	sed -e 's|.*ckbs-\([0-9.]\{10\}\).*|\1|'`
echo "doc_version = $doc_version"
#
if [ "$doc_version" != "$stable_version" ]
then
	echo "new_stable.sh: run ./copy_doc.sh to bring doc_version up to date"
	echo "new_stable.sh: stable_version = $stable_version"
	echo "new_stable.sh: doc_version    = $doc_version"
	exit 1
fi
#
# web address for trunk, stable, release, and documentation
repository="https://projects.coin-or.org/svn/CoinBazaar"
rep_trunk="$repository/projects/ckbs/trunk"
rep_stable="$repository/projects/ckbs/stable/$stable_version"
rep_release="$repository/projects/ckbs/releases/$stable_version.0"
rep_html="$repository/html/ckbs"
#
# create the new stable version
msg="copy ckbs/trunk at revision $rev_trunk to ckbs/stable/$stable_version"
echo "svn copy $rep_trunk $rep_stable -m \"$msg\""
svn copy  $rep_trunk $rep_stable -m "$msg"
#
# Add documentation to the stable version
msg="[ckbs/stable] Add documentation corresponding to this stable version"
echo "svn copy $rep_html $rep_stable/doc -m \"$msg\""
svn copy $rep_html $rep_stable/doc -m "$msg"
#
# create the new release version
msg="copy ckbs/stable/$stable_version to ckbs/releases/$stable_version.0"
svn copy $rep_stable $rep_release -m "$msg"
#
echo ""
echo "To check out this stable version use the command:"
echo "	svn checkout \\"
echo "$rep_stable"
echo "To check out this release version use the command:"
echo "	svn checkout \\"
echo "$rep_release"
