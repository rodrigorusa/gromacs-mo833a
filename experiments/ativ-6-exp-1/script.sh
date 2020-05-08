# Compiling
echo "Compiling..."
cd ../../
[ ! -d "build/" ] && mkdir -p "build/"
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=OFF -DBUILD_SHARED_LIBS=OFF
make -j8

echo "Creating simulation..."
cd ../experiments/ativ-6-exp-1

# Criar a topologia do ambiente que irá ser simulado:
echo "15" | ./../../build/bin/gmx pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build/bin/gmx editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build/bin/gmx solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build/bin/gmx grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
echo "13" | ./../../build/bin/gmx genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build/bin/gmx grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr

echo "Running..."
./../../build/bin/gmx mdrun -v -deffnm em >> paramount_time.csv

find . -maxdepth 1 -not -type d ! -name "6LVN.pdb" ! -name "ions.mdp" ! -name "script.sh" ! -name "README.md" ! -name "paramount_time.csv" -type f -exec rm -f {} +
