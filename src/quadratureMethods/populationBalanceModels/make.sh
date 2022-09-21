echo "start making..."
wmake -j
echo "moving"
rm /home/ubuntu/OpenFOAM/ubuntu-v2206/platforms/linux64GccDPInt32Opt/lib/libpopulationBalance.so
cp /root/OpenFOAM/root-v2206/platforms/linux64GccDPInt32Opt/lib/libpopulationBalance.so /home/ubuntu/OpenFOAM/ubuntu-v2206/platforms/linux64GccDPInt32Opt/lib/

