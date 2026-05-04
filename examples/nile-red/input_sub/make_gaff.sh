antechamber -i nile-red.sdf -fi sdf -o nile-red_gaff2.mol2 -fo mol2 -c bcc -nc 0 -rn QUC -at gaff2

parmchk2 -i nile-red_gaff2.mol2 -f mol2 -o nile-red_gaff2.frcmod -s gaff2

antechamber -i nile-red_gaff2.mol2 -fi mol2 -o nile-red.pdb -fo pdb
