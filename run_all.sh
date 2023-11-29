make clean
make -j8
./Fake
./Fake_lxy

root -l -b -q fitks_diffun_John.cc >&log_pt_fits&
root -l -b -q fitks_lxy.cc >&log_lxy_fits&
