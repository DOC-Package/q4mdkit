#!/bin/sh

antechamber -i picene_template.mol2 -fi mol2 -o picene_gaff2.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn PIC
parmchk2 -i picene_gaff2.mol2 -f mol2 -o picene.frcmod -s gaff2
tleap -f tleap.in
#cpptraj -i setbox.cpptraj