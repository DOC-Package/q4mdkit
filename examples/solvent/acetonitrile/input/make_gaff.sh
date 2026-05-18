antechamber -i acetonitrile.sdf -fi sdf -o acetonitrile_gaff2.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn ACE

parmchk2 -i acetonitrile_gaff2.mol2 -f mol2 -o acetonitrile_gaff2.frcmod -s gaff2

antechamber -i acetonitrile_gaff2.mol2 -fi mol2 -o acetonitrile.pdb -fo pdb

#tleap -f leap_acetonitrile_box.in