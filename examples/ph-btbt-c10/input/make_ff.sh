#!/bin/sh

antechamber -i ph-btbt-c10_template.mol2 -fi mol2 -o ph-btbt-c10_gaff2.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn BTB
parmchk2 -i ph-btbt-c10_gaff2.mol2 -f mol2 -o ph-btbt-c10.frcmod -s gaff2
tleap -f tleap.in
#cpptraj -i setbox.cpptraj