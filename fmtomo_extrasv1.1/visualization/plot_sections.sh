#!/bin/bash

cd gmt_files
path="$(cat data_path)"
while read -r spacing; do echo ""; done < Spacing
n=1
while read line; do
coor[n]=$line
n=$((n+1))
done < Area_boundaries
Xmin=${coor[1]}
Xmax=${coor[2]}
Ymin=${coor[3]}
Ymax=${coor[4]}
cp ../gmtslicei_empty.in $path/invert_p/gmtplot/gmtsliceiP.in
cp ../gmtslicei_empty.in $path/invert_s/gmtplot/gmtsliceiS.in
cp * $path/invert_p/gmtplot/
cp * $path/invert_s/gmtplot/

echo "Starting plot session. It needs a lot of time. Please be patient."
echo "On the terminal it will appear the actual plot."
echo ""

# Plot P model
echo "Starting with P model"
cd $path/invert_p/gmtplot/
mkdir -p P_images
mkdir -p P_images/absolute
mkdir -p P_images/relative
mv gmtslicei.in gmtslicei.in.or
# absolute velocity
velocity=velabsolute.cpt
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/1/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NZ;
do
    echo "Absolute - Horizontal slice $NZ"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    gmtslicei
    #plotd
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < bounddp.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj=$(cat Project)
    gmt begin HOR pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs8c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvd.z -Ggrid2dvd.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray
            gmt coast $bounds $proj -B -Na0.5 -B+glightblue -W0.3p,black
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity -t50
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt subplot set # second subplot: rays
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray
            gmt coast $bounds $proj -B -Na0.5 -B+glightblue -W0.3p,black
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity -t50
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt psxy raysd.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesd.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversd.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv HOR.ps P_images/absolute/HOR_$NZ.ps
    mv HOR.pdf P_images/absolute/HOR_$NZ.pdf
    mv HOR.png P_images/absolute/HOR_$NZ.png
    rm gmtslicei.in
done < DEPTH
rm gmtsliceiZ.in
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/1/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NY;
do
    echo "Absolute - EW slice $NY"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LATY|$NY|g" gmtslicei.in
    gmtslicei
    #plotew
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundns.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin EW pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvns.z -Ggrid2dvns.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xmin/$NY/$Xmax/$NY -N | \
            awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xmin/$NY -E$Xmax/$NY -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysns.xy $bounds $proj -W1p,black
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt psxy sourcesns.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversns.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv EW.ps P_images/absolute/EW_$NY.ps
    mv EW.pdf P_images/absolute/EW_$NY.pdf
    mv EW.png P_images/absolute/EW_$NY.png
    rm gmtslicei.in
done < LATY
rm gmtsliceiZ.in
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/1/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
while read NX;
do
    echo "Absolute - SN slice $NX"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LONGY|$NX|g" gmtslicei.in
    gmtslicei
    #plotsn
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundew.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin SN pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvew.z -Ggrid2dvew.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$NX/$Ymin/$NX/$Ymax -N | \
            awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$NX/$Ymin -E$NX/$Ymax -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysew.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesew.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversew.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv SN.ps P_images/absolute/SN_$NX.ps
    mv SN.pdf P_images/absolute/SN_$NX.pdf
    mv SN.png P_images/absolute/SN_$NX.png
    rm gmtslicei.in
done < LONGY
rm gmtsliceiZ.in
while read -r Ystart Xstart Yend Xend SezLen SezName;
do
    echo "Absolute - Vertical slice $SezName"
    cp gmtsliceiP.in gmtslicei.in
    sed -i 's/VEL/0/g' gmtslicei.in
    sed -i 's/GCS/1/g' gmtslicei.in
    sed -i 's/HOR/0/g' gmtslicei.in
    sed -i 's/EWS/0/g' gmtslicei.in
    sed -i 's/NSS/0/g' gmtslicei.in
    sed -i "s|STARTLAT|$Ystart|g" gmtslicei.in
    sed -i "s|STARTLON|$Xstart|g" gmtslicei.in
    sed -i "s|ENDLAT|$Yend|g" gmtslicei.in
    sed -i "s|ENDLON|$Xend|g" gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    sed -i "s|LATY|$Ymin|g" gmtslicei.in
    sed -i "s|LONGY|$Xmin|g" gmtslicei.in
    gmtslicei
    #plotgc
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundgc.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin GC pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvgc.z -Ggrid2dvgc.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xstart/$Ystart/$Xend/$Yend -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xstart/$Ystart -E$Xend/$Yend -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysgc.xy $bounds $proj -W1p,black
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt psxy sourcesgc.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversgc.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv GC.ps P_images/absolute/GC_$SezName.ps
    mv GC.pdf P_images/absolute/GC_$SezName.pdf
    mv GC.png P_images/absolute/GC_$SezName.png
    rm gmtslicei.in
done < Section_coordinates
# relative velocity
velocity=velrelative.cpt
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/1/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NZ;
do
    echo "Relative - Horizontal slice $NZ"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    gmtslicei
    #plotd
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < bounddp.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj=$(cat Project)
    gmt begin HOR pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs8c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvd.z -Ggrid2dvd.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt psxy raysd.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesd.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversd.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
        gmt subplot end
    gmt end
    mv HOR.ps P_images/absolute/HOR_$NZ.ps
    mv HOR.pdf P_images/absolute/HOR_$NZ.pdf
    mv HOR.png P_images/absolute/HOR_$NZ.png
    rm gmtslicei.in
done < DEPTH
rm gmtsliceiZ.in
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/1/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NY;
do
    echo "Relative - EW slice $NY"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LATY|$NY|g" gmtslicei.in
    gmtslicei
    #plotew
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundns.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin EW pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvns.z -Ggrid2dvns.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xmin/$NY/$Xmax/$NY -N | \
            awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xmin/$NY -E$Xmax/$NY -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysns.xy $bounds $proj -W1p,black
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt psxy sourcesns.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversns.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv EW.ps P_images/absolute/EW_$NY.ps
    mv EW.pdf P_images/absolute/EW_$NY.pdf
    mv EW.png P_images/absolute/EW_$NY.png
    rm gmtslicei.in
done < LATY
rm gmtsliceiZ.in
cp gmtsliceiP.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/1/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
while read NX;
do
    echo "Relative - SN slice $NX"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LONGY|$NX|g" gmtslicei.in
    gmtslicei
    #plotsn
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundew.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin SN pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvew.z -Ggrid2dvew.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$NX/$Ymin/$NX/$Ymax -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$NX/$Ymin -E$NX/$Ymax -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysew.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesew.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversew.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv SN.ps P_images/absolute/SN_$NX.ps
    mv SN.pdf P_images/absolute/SN_$NX.pdf
    mv SN.png P_images/absolute/SN_$NX.png
    rm gmtslicei.in
done < LONGY
rm gmtsliceiZ.in
while read -r Ystart Xstart Yend Xend SezLen SezName;
do
    echo "Relative - Vertical slice $SezName"
    cp gmtsliceiP.in gmtslicei.in
    sed -i 's/VEL/0/g' gmtslicei.in
    sed -i 's/GCS/1/g' gmtslicei.in
    sed -i 's/HOR/0/g' gmtslicei.in
    sed -i 's/EWS/0/g' gmtslicei.in
    sed -i 's/NSS/0/g' gmtslicei.in
    sed -i "s|STARTLAT|$Ystart|g" gmtslicei.in
    sed -i "s|STARTLON|$Xstart|g" gmtslicei.in
    sed -i "s|ENDLAT|$Yend|g" gmtslicei.in
    sed -i "s|ENDLON|$Xend|g" gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    sed -i "s|LATY|$Ymin|g" gmtslicei.in
    sed -i "s|LONGY|$Xmin|g" gmtslicei.in
    gmtslicei
    #plotgc
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundgc.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin GC pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvgc.z -Ggrid2dvgc.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xstart/$Ystart/$Xend/$Yend -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xstart/$Ystart -E$Xend/$Yend -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysgc.xy $bounds $proj -W1p,black
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt psxy sourcesgc.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversgc.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv GC.ps P_images/absolute/GC_$SezName.ps
    mv GC.pdf P_images/absolute/GC_$SezName.pdf
    mv GC.png P_images/absolute/GC_$SezName.png
    rm gmtslicei.in
done < Section_coordinates
mv gmtslicei.in.or gmtslicei.in

# Plot S model
echo "Changing to S model"
cd ../../invert_s/gmtplot/
mkdir -p S_model
mkdir -p S_model/absolute
mkdir -p S_model/relative
mv gmtslicei.in gmtslicei.in.or
# absolute velocity
velocity=velabsolute.cpt
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/1/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NZ;
do
    echo "Absolute - Horizontal slice $NZ"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    gmtslicei
    #plotd
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < bounddp.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj=$(cat Project)
    gmt begin HOR pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs8c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvd.z -Ggrid2dvd.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt psxy raysd.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesd.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversd.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
        gmt subplot end
    gmt end
    mv HOR.ps S_images/absolute/HOR_$NZ.ps
    mv HOR.pdf S_images/absolute/HOR_$NZ.pdf
    mv HOR.png S_images/absolute/HOR_$NZ.png
    rm gmtslicei.in
done < DEPTH
rm gmtsliceiZ.in
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/1/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NY;
do
    echo "Absolute - EW slice $NY"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LATY|$NY|g" gmtslicei.in
    gmtslicei
    #plotew
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundns.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin EW pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvns.z -Ggrid2dvns.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xmin/$NY/$Xmax/$NY -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xmin/$NY -E$Xmax/$NY -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysns.xy $bounds $proj -W1p,black
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt psxy sourcesns.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversns.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv EW.ps S_images/absolute/EW_$NY.ps
    mv EW.pdf S_images/absolute/EW_$NY.pdf
    mv EW.png S_images/absolute/EW_$NY.png
    rm gmtslicei.in
done < LATY
rm gmtsliceiZ.in
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/1/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
while read NX;
do
    echo "Absolute - SN slice $NX"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LONGY|$NX|g" gmtslicei.in
    gmtslicei
    #plotsn
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundew.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin SN pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvew.z -Ggrid2dvew.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$NX/$Ymin/$NX/$Ymax -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$NX/$Ymin -E$NX/$Ymax -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysew.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesew.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversew.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv SN.ps S_images/absolute/SN_$NX.ps
    mv SN.pdf S_images/absolute/SN_$NX.pdf
    mv SN.png S_images/absolute/SN_$NX.png
    rm gmtslicei.in
done < LONGY
rm gmtsliceiZ.in
while read -r Ystart Xstart Yend Xend SezLen SezName;
do
    echo "Absolute - Vertical slice $SezName"
    cp gmtsliceiS.in gmtslicei.in
    sed -i 's/VEL/0/g' gmtslicei.in
    sed -i 's/GCS/1/g' gmtslicei.in
    sed -i 's/HOR/0/g' gmtslicei.in
    sed -i 's/EWS/0/g' gmtslicei.in
    sed -i 's/NSS/0/g' gmtslicei.in
    sed -i "s|STARTLAT|$Ystart|g" gmtslicei.in
    sed -i "s|STARTLON|$Xstart|g" gmtslicei.in
    sed -i "s|ENDLAT|$Yend|g" gmtslicei.in
    sed -i "s|ENDLON|$Xend|g" gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    sed -i "s|LATY|$Ymin|g" gmtslicei.in
    sed -i "s|LONGY|$Xmin|g" gmtslicei.in
    gmtslicei
    #plotgc
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundgc.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin GC pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvgc.z -Ggrid2dvgc.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xstart/$Ystart/$Xend/$Yend -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xstart/$Ystart -E$Xend/$Yend -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysgc.xy $bounds $proj -W1p,black
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt psxy sourcesgc.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversgc.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv GC.ps S_images/absolute/GC_$SezName.ps
    mv GC.pdf S_images/absolute/GC_$SezName.pdf
    mv GC.png S_images/absolute/GC_$SezName.png
    rm gmtslicei.in
done < Section_coordinates
# relative velocity
velocity=velrelative.cpt
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/1/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NZ;
do
    echo "Relative - Horizontal slice $NZ"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    gmtslicei
    #plotd
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < bounddp.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj=$(cat Project)
    gmt begin HOR pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs8c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvd.z -Ggrid2dvd.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvd.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx8.5c/1.5c+w5c/0.3c
            gmt grdimage @earth_gebco_02m $bounds $proj -Cgray -t80
            gmt psxy raysd.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesd.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversd.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
            gmt coast $bounds $proj -B -Glightgray -Na0.5 -B+glightblue -W0.3p,black
        gmt subplot end
    gmt end
    mv HOR.ps S_images/absolute/HOR_$NZ.ps
    mv HOR.pdf S_images/absolute/HOR_$NZ.pdf
    mv HOR.png S_images/absolute/HOR_$NZ.png
    rm gmtslicei.in
done < DEPTH
rm gmtsliceiZ.in
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/1/g' gmtsliceiZ.in
sed -i 's/NSS/0/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LONGY|$Xmin|g" gmtsliceiZ.in
while read NY;
do
    echo "Relative - EW slice $NY"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LATY|$NY|g" gmtslicei.in
    gmtslicei
    #plotew
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundns.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin EW pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvns.z -Ggrid2dvns.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xmin/$NY/$Xmax/$NY -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xmin/$NY -E$Xmax/$NY -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvns.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysns.xy $bounds $proj -W1p,black
            gmt psxy intns.xy $bounds $proj -W1p,black
            gmt psxy sourcesns.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversns.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv EW.ps S_images/absolute/EW_$NY.ps
    mv EW.pdf S_images/absolute/EW_$NY.pdf
    mv EW.png S_images/absolute/EW_$NY.png
    rm gmtslicei.in
done < LATY
rm gmtsliceiZ.in
cp gmtsliceiS.in gmtsliceiZ.in
sed -i 's/VEL/0/g' gmtsliceiZ.in
sed -i 's/HOR/0/g' gmtsliceiZ.in
sed -i 's/EWS/0/g' gmtsliceiZ.in
sed -i 's/NSS/1/g' gmtsliceiZ.in
sed -i 's/GCS/0/g' gmtsliceiZ.in
sed -i "s|STARTLAT|$Ymin|g" gmtsliceiZ.in
sed -i "s|STARTLON|$Xmin|g" gmtsliceiZ.in
sed -i "s|ENDLAT|$Ymax|g" gmtsliceiZ.in
sed -i "s|ENDLON|$Xmax|g" gmtsliceiZ.in
sed -i "s|DEPTH|$NZ|g" gmtsliceiZ.in
sed -i "s|LATY|$Ymin|g" gmtsliceiZ.in
while read NX;
do
    echo "Relative - SN slice $NX"
    cp gmtsliceiZ.in gmtslicei.in
    sed -i "s|LONGY|$NX|g" gmtslicei.in
    gmtslicei
    #plotsn
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundew.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin SN pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvew.z -Ggrid2dvew.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$NX/$Ymin/$NX/$Ymax -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$NX/$Ymin -E$NX/$Ymax -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvew.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysew.xy $bounds $proj -W1p,black
            gmt psxy intew.xy $bounds $proj -W1p,black
            gmt psxy sourcesew.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversew.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv SN.ps S_images/absolute/SN_$NX.ps
    mv SN.pdf S_images/absolute/SN_$NX.pdf
    mv SN.png S_images/absolute/SN_$NX.png
    rm gmtslicei.in
done < LONGY
rm gmtsliceiZ.in
while read -r Ystart Xstart Yend Xend SezLen SezName;
do
    echo "Relative - Vertical slice $SezName"
    cp gmtsliceiS.in gmtslicei.in
    sed -i 's/VEL/0/g' gmtslicei.in
    sed -i 's/GCS/1/g' gmtslicei.in
    sed -i 's/HOR/0/g' gmtslicei.in
    sed -i 's/EWS/0/g' gmtslicei.in
    sed -i 's/NSS/0/g' gmtslicei.in
    sed -i "s|STARTLAT|$Ystart|g" gmtslicei.in
    sed -i "s|STARTLON|$Xstart|g" gmtslicei.in
    sed -i "s|ENDLAT|$Yend|g" gmtslicei.in
    sed -i "s|ENDLON|$Xend|g" gmtslicei.in
    sed -i "s|DEPTH|$NZ|g" gmtslicei.in
    sed -i "s|LATY|$Ymin|g" gmtslicei.in
    sed -i "s|LONGY|$Xmin|g" gmtslicei.in
    gmtslicei
    #plotgc
    n=1
    while read line; do
    bds[n]=$line
    n=$((n+1))
    done < boundgc.gmt
    bounds="-R${bds[1]}/${bds[2]}/${bds[3]}/${bds[4]}"
    proj="-JX17c/8c"
    gmt begin GC pdf,png,ps
        gmt set FONT_ANNOT 10p,Helvetica
        gmt subplot begin 2x1 -Fs18c/8c
            gmt subplot set # first plot: velocity only
            gmt xyz2grd grid2dvgc.z -Ggrid2dvgc.grd -I${bds[5]}+/${bds[6]}+ $bounds -ZLB
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt grdtrack -G@earth_relief_01m -E$Xstart/$Ystart/$Xend/$Yend -N | awk '{printf "%s\t%.3f\n", $0, $3/1000}' > topography
            gmt project topography -C$Xstart/$Ystart -E$Xend/$Yend -Q -W-30/30 -i0,1,3 > topo_track
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt subplot set # second subplot: rays
            gmt grdimage grid2dvgc.grd $bounds $proj -Bxa1.0f0.5 -Bya40f10 -C$velocity
            gmt psscale -C$velocity -Ba0.5f0.1 -Dx17.5c/1.5c+w5c/0.3c
            gmt plot topo_track $bounds $proj -W2p,black -i3,2
            gmt psxy raysgc.xy $bounds $proj -W1p,black
            gmt psxy intgc.xy $bounds $proj -W1p,black
            gmt psxy sourcesgc.xy $bounds $proj -Si0.40c -Gblue -N -W1p,black
            gmt psxy receiversgc.xy $bounds $proj -Sa0.20c -Gred -N -W0.3p,black
        gmt subplot end
    gmt end
    mv GC.ps S_images/absolute/GC_$SezName.ps
    mv GC.pdf S_images/absolute/GC_$SezName.pdf
    mv GC.png S_images/absolute/GC_$SezName.png
    rm gmtslicei.in
done < Section_coordinates
mv gmtslicei.in.or gmtslicei.in
