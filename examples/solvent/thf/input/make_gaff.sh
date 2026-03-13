antechamber -i thf.sdf -fi sdf -o thf_gaff2.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn THF

parmchk2 -i thf_gaff2.mol2 -f mol2 -o thf_gaff2.frcmod -s gaff2

antechamber -i thf_gaff2.mol2 -fi mol2 -o thf.pdb -fo pdb

#tleap -f leap_thf_box.in