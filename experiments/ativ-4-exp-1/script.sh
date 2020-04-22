# Compiling Release mode with OpenMP, MPI and GPU
echo "Compiling release build with OpenMP..."
cd ../../
[ ! -d "build-release-openmp/" ] && mkdir -p "build-release-openmp/"
cd build-release-openmp
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=/usr/bin/gcc-6 -DCMAKE_CXX_COMPILER=/usr/bin/g++-6
make -j8
make check -j8

# Compiling Release mode without OpenMP, MPI and GPU
echo "Compiling release build without OpenMP..."
cd ../
[ ! -d "build-release-no-openmp/" ] && mkdir -p "build-release-no-openmp/"
cd build-release-no-openmp
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_THREAD_MPI=OFF -DGMX_MPI=OFF -DGMX_OPENMP=OFF -DGMX_GPU=OFF -DCMAKE_C_COMPILER=/usr/bin/gcc-6 -DCMAKE_CXX_COMPILER=/usr/bin/g++-6
make -j8
make check -j8

echo "Creating simulation with OpenMP..."
cd ../experiments/ativ-4-exp-1

# Criar a topologia do ambiente que irá ser simulado:
./../../build-release-openmp/bin/gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build-release-openmp/bin/gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build-release-openmp/bin/gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build-release-openmp/bin/gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
./../../build-release-openmp/bin/gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build-release-openmp/bin/gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr

echo "Profiling with OpenMP..."
/usr/bin/time -o time_exec.txt -p sudo perf record ./../../build-release-openmp/bin/gmx mdrun -v -deffnm em

find . -maxdepth 1 -not -type d ! -name "6LVN.pdb" ! -name "ions.mdp" ! -name "script.sh" ! -name "perf.data" ! -name "time_exec.txt" -type f -exec rm -f {} +

[ ! -d "perf-openmp/" ] && mkdir -p "perf-openmp/"
sudo perf report > perf_report.txt
sudo mv perf_report.txt perf-openmp/perf_report.txt
sudo mv perf.data perf-openmp/perf.data
mv time_exec.txt perf-openmp/time_exec.txt

echo "Creating simulation without OpenMP..."

# Criar a topologia do ambiente que irá ser simulado:
./../../build-release-no-openmp/bin/gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build-release-no-openmp/bin/gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build-release-no-openmp/bin/gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build-release-no-openmp/bin/gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
./../../build-release-no-openmp/bin/gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build-release-no-openmp/bin/gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr

echo "Profiling without OpenMP..."
/usr/bin/time -o time_exec.txt -p sudo perf record ./../../build-release-no-openmp/bin/gmx mdrun -v -deffnm em

find . -maxdepth 1 -not -type d ! -name "6LVN.pdb" ! -name "ions.mdp" ! -name "script.sh" ! -name "perf.data" ! -name "time_exec.txt" ! -name "perf_report.txt" -type f -exec rm -f {} +

[ ! -d "perf-no-openmp/" ] && mkdir -p "perf-no-openmp/"
sudo perf report > perf_report.txt
sudo mv perf_report.txt perf-no-openmp/perf_report.txt
sudo mv perf.data perf-no-openmp/perf.data
mv time_exec.txt perf-no-openmp/time_exec.txt
