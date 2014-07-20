#! /bin/sh
# -------------------------------------------------------------------
# ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
# Authors: Bradlely Bell:        bradbell at washington dot edu
#          Gianluigi Pillonetto: giapi at dei dot unipd dot it
# License: GNU General Public License Version 2
# -------------------------------------------------------------------
#
# exit on error
set -e
#
for dir in example test
do
	cd $dir
	cat << EOF > test.$$
more off
if all_ok
	exit(0)
else
	exit(1)
end
EOF
	echo "test.sh: running $dir/all_ok.m"
	if ! octave --silent test.$$
	then
		rm test.$$
		echo "test.sh: $dir/all_ok failed"
		exit 1
	fi
	rm test.$$
	cd ..
done
exit 0
