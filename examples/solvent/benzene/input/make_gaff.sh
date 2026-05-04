antechamber -i benzene.sdf -fi sdf -o benzene_gaff2.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn DMS

parmchk2 -i benzene_gaff2.mol2 -f mol2 -o benzene_gaff2.frcmod -s gaff2

antechamber -i benzene_gaff2.mol2 -fi mol2 -o benzene.pdb -fo pdb

#tleap -f leap_benzene_box.in