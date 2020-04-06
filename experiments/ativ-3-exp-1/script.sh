# Compiling Release mode
echo "Compiling release build..."
cd ../../
[ ! -d "build-release/" ] && mkdir -p "build-release/"
cd build-release
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=/usr/bin/gcc-6 -DCMAKE_CXX_COMPILER=/usr/bin/g++-6
make -j8
make check -j8

# Compiling Debug mode
echo "Compiling debug build..."
cd ../
[ ! -d "build-debug/" ] && mkdir -p "build-debug/"
cd build-debug
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=/usr/bin/gcc-6 -DCMAKE_CXX_COMPILER=/usr/bin/g++-6 -DCMAKE_BUILD_TYPE=Debug
make -j8
make check -j8

echo "Creating simulation..."
cd ../experiments/ativ-3-exp-1/

# Criar a topologia do ambiente que irá ser simulado:
./../../build-release/bin/gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build-release/bin/gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build-release/bin/gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build-release/bin/gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
./../../build-release/bin/gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build-release/bin/gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr

# Executar simulação:
echo "Executing release build..."
for i in $(seq 1 99);
do
	./../../build-release/bin/gmx mdrun -v -deffnm em | grep "[MO833]" >> time_exec_release.csv
done

find . ! -name "6LVN.pdb" ! -name "ions.mdp" ! -name "script.sh" ! -name "time_exec_release.csv" -type f -exec rm -f {} +

# Criar a topologia do ambiente que irá ser simulado:
./../../build-debug/bin/gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build-debug/bin/gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build-debug/bin/gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build-debug/bin/gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
./../../build-debug/bin/gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build-debug/bin/gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr

echo "Executing debug build..."
for i in $(seq 1 99);
do
	./../../build-debug/bin/gmx mdrun -v -deffnm em | grep "[MO833]" >> time_exec_debug.csv
done

find . ! -name "6LVN.pdb" ! -name "ions.mdp" ! -name "script.sh" ! -name "time_exec_release.csv" ! -name "time_exec_debug.csv" -type f -exec rm -f {} +
