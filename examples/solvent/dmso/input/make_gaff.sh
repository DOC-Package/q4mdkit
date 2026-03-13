antechamber -i dmso.sdf -fi sdf -o dmso_gaff2.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn DMS

parmchk2 -i dmso_gaff2.mol2 -f mol2 -o dmso_gaff2.frcmod -s gaff2

antechamber -i dmso_gaff2.mol2 -fi mol2 -o dmso.pdb -fo pdb

#tleap -f leap_dmso_box.in