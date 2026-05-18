antechamber -i pdi.sdf -fi sdf -o pdi_gaff2.mol2 -fo mol2 -c bcc -nc 0 -rn PDI -at gaff2

parmchk2 -i pdi_gaff2.mol2 -f mol2 -o pdi_gaff2.frcmod -s gaff2

antechamber -i pdi_gaff2.mol2 -fi mol2 -o pdi.pdb -fo pdb
