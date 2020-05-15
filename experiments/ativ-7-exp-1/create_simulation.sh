echo "Creating simulation..."

# Criar a topologia do ambiente que irá ser simulado:
echo "15" | ./../../build/bin/gmx_mpi pdb2gmx -f 6LVN.pdb -o 6LVN_processed.gro -water spce

# Definir a “caixa” no qual a molécula, os ions e a água irá estar:
./../../build/bin/gmx_mpi editconf -f 6LVN_processed.gro -o 6LVN_newbox.gro -c -d 1.0 -bt cubic

# Adicionar o solvente (água) na caixa:
./../../build/bin/gmx_mpi solvate -cp 6LVN_newbox.gro -cs spc216.gro -o 6LVN_solv.gro -p topol.top

# Adicionar os ions na caixa:
./../../build/bin/gmx_mpi grompp -f ions.mdp -c 6LVN_solv.gro -p topol.top -o ions.tpr
echo "13" | ./../../build/bin/gmx_mpi genion -s ions.tpr -o 6LVN_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Gerar a simulação:
./../../build/bin/gmx_mpi grompp -f ions.mdp -c 6LVN_solv_ions.gro -p topol.top -o em.tpr