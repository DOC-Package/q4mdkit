#!/bin/sh

antechamber -i spiro-ometad_template.mol2 -fi mol2 -o spiro-ometad_gaff2.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn SPO
parmchk2 -i spiro-ometad_gaff2.mol2 -f mol2 -o spiro-ometad.frcmod -s gaff2
tleap -f tleap.in
#cpptraj -i setbox.cpptraj