#! /bin/bash
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# ----------------------------------------------------------------------------
stable_version="0.20130204"
release_number="1"
msg="[ckbs/releases] add documentation (missing from 0.20100228.0 release)."
# -------------------------------------------------------------------------
repository="https://projects.coin-or.org/svn/CoinBazaar/projects/ckbs"
rep_stable=$repository/stable/$stable_version
rep_release=$repository/releases/$stable_version.$release_number
#
echo "svn $rep_stable $rep_release -m \"$msg\""
if ! svn copy $rep_stable $rep_release -m "$msg"
then
	echo "new_stable.sh: svn copy failed"
fi
exit 0
