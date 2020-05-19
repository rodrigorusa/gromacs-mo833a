# Compiling
echo "Compiling..."
cd ../../
[ ! -d "build/" ] && mkdir -p "build/"
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_SHARED_LIBS=OFF
make -j4