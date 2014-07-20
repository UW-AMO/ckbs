#! /bin/sh
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# -------------------------------------------------------------------
#
svn status | sed \
	-e '/^[?] *commit.sh$/d' \
	-e '/^[?] *doc$/d' \
	-e '/^[?] *junk$/d' \
	-e '/^[?] *junk.sh$/d' \
	-e '/^[?] *omhelp.log$/d' \
	-e '/^[?] *omhelp-[0-9][0-9]-[0-9][0-9]-[0-9][0-9]$/d' \
	-e '/^[?] *OMhelp.unix.tar.gz$/d' \
	-e '/^[?] *temp$/d'
	
