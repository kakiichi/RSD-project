# ! /usr/bin/tcsh

# REF
#./src/IGMtransmission.exe -in input/param_01.in -out output/IGMtransmission_01.dat | tee log/IGMtransmission_01.log
#./src/IGMtransmission.exe -in input/param_01.in -out output/IGMtransmission_01_100kpc.dat | tee log/IGMtransmission_01.log
#./src/IGMtransmission.exe -in input/param_01.in -out output/IGMtransmission_01_30kpc.dat | tee log/IGMtransmission_01.log

# DENS
#./src/RSD.exe -in input/param_02.in -out output/RSD_02.dat > log/RSD_02.log &
#./src/RSD.exe -in input/param_03.in -out output/RSD_03.dat > log/RSD_03.log &

# LyC
#./src/RSD.exe -in input/param_04.in -out output/RSD_04.dat > log/RSD_04.log &
#./src/RSD.exe -in input/param_05.in -out output/RSD_05.dat > log/RSD_05.log &

# inflow
#./src/RSD.exe -in input/param_06.in -out output/RSD_06.dat > log/RSD_06.log &
#./src/RSD.exe -in input/param_07.in -out output/RSD_07.dat > log/RSD_07.log &

# outflow
#./src/RSD.exe -in input/param_08.in -out output/RSD_08.dat > log/RSD_08.log &
#./src/RSD.exe -in input/param_09.in -out output/RSD_09.dat > log/RSD_09.log &
#./src/RSD.exe -in input/param_10.in -out output/RSD_10.dat > log/RSD_10.log &
./src/IGMtransmission.exe -in input/param_08.in -out output/IGMtransmission_08_30kpc.dat | tee log/IGMtransmission_08.log
./src/IGMtransmission.exe -in input/param_09.in -out output/IGMtransmission_09_30kpc.dat | tee log/IGMtransmission_09.log
./src/IGMtransmission.exe -in input/param_10.in -out output/IGMtransmission_10_30kpc.dat | tee log/IGMtransmission_10.log

# velocity dispersion
#./src/IGMtransmission.exe -in input/param_11.in -out output/IGMtransmission_11_30kpc.dat | tee log/IGMtransmission_11.log
#./src/IGMtransmission.exe -in input/param_12.in -out output/IGMtransmission_12_30kpc.dat | tee log/IGMtransmission_12.log


