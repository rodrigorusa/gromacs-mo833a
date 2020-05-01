# Compiling
echo "Compiling..."
cd ../../
[ ! -d "build/" ] && mkdir -p "build/"
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_MPI=ON
make