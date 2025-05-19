#!/bin/bash

echo "preparing files for tomography"
rm frechet.dat interfaces* mtimes.dat otimes* rays.dat sandr.dat sources* src.dat vgrid* vpvsresids.out
ksh cp_script
echo "performing the Vp/Vs tomography"
tomo3dps
