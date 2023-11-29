make clean
make -j8
./Fake

root -l -b -q fitks.cc >&log_fits&
