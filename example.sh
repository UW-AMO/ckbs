#! /bin/sh
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# -------------------------------------------------------------------
#
cd example
cat << EOF > example.$$
if all_ok
	exit(0)
else
	exit(1)
end
EOF
if ! octave --silent example.$$
then
	rm example.$$
	echo "example.sh: example/all_ok failed"
	exit 1
fi
rm example.$$
exit 0
