antechamber -i acetone.sdf -fi sdf -o acetone_gaff2.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn ACE

parmchk2 -i acetone_gaff2.mol2 -f mol2 -o acetone_gaff2.frcmod -s gaff2

antechamber -i acetone_gaff2.mol2 -fi mol2 -o acetone.pdb -fo pdb

#tleap -f leap_acetone_box.in