#!/bin/sh

antechamber -i pentacene_template.mol2 -fi mol2 -o pentacene_gaff2.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn PEN
parmchk2 -i pentacene_gaff2.mol2 -f mol2 -o pentacene.frcmod -s gaff2
tleap -f tleap.in
cpptraj -i setbox.cpptraj