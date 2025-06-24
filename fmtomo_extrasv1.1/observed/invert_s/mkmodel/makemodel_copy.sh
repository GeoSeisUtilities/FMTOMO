#!/bin/bash

echo "preparing earthquake data"
rm sources.in sourcesref.in receivers.in receiversref.in otimes.dat otimesref.dat vgrids.in vgridsref.in interfaces.in interfacesref.in src.dat
obsdatai
cp sources.in sourcesref.in
cp receivers.in receiversref.in
cp otimes.dat otimesref.dat
echo "preparing model geometry"
grid3dgi
cp vgrids.in vgridsref.in
cp vgridsref.in reference_velocity.txt
cp interfaces.in interfacesref.in
echo "copying files"
cp sources* ../
cp vgrids* ../
cp interfaces* ../
cp otimes* ../
cp propgrid.in ../
cp receivers* ../
cp src.dat ../
