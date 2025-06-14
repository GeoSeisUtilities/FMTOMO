#!/bin/bash

echo "Removing old files"
rm arrivals.dat arrtimes.dat frechet.dat gridsave.in
echo "Creating gridsave.in from $1"
n=$(head -n 1 $1)
touch gridsave.in
m=1
for i in $(seq 1 $n)
do
    echo "$m 1" >> gridsave.in
    echo "1" >> gridsave.in
    echo "1" >> gridsave.in
    m=$(($m + 1))
done
echo "Calulating arrivals"
fm3d
echo "Setting synthetic model"
cd mkmodel
cp ../arrivals.dat ./
synthdatai
cp src.dat ../
cp otimes.dat otimesref.dat
cp otimesref.dat ../
