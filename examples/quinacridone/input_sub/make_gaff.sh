antechamber -i quinacridone.sdf -fi sdf -o quinacridone_gaff2.mol2 -fo mol2 -c bcc -nc 0 -rn QUC -at gaff2

parmchk2 -i quinacridone_gaff2.mol2 -f mol2 -o quinacridone_gaff2.frcmod -s gaff2

antechamber -i quinacridone_gaff2.mol2 -fi mol2 -o quinacridone.pdb -fo pdb
