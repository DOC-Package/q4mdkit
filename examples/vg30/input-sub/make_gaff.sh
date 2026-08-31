#!/bin/sh
set -eu

xyz_file="vg30-6-cl.xyz"
gcrt_file=$(mktemp "${TMPDIR:-/tmp}/vg30-6-cl.XXXXXX")
trap 'rm -f "$gcrt_file"' EXIT HUP INT TERM

awk '
    NR == 1 {
        natoms = $1 + 0
        if (natoms <= 0) {
            print "Invalid XYZ atom count" > "/dev/stderr"
            exit 1
        }
        print "%chk=vg30-6-cl.chk"
        print "# AM1"
        print ""
        print "VG30-6-Cl cation"
        print ""
        print "1 1"
        next
    }
    NR == 2 { next }
    atom_count < natoms {
        if (NF < 4) {
            print "Invalid XYZ coordinate line " NR > "/dev/stderr"
            exit 1
        }
        printf "%-2s %16.8f %16.8f %16.8f\n", $1, $2, $3, $4
        atom_count++
    }
    END {
        if (atom_count != natoms) {
            print "XYZ atom count does not match coordinate lines" > "/dev/stderr"
            exit 1
        }
        print ""
    }
' "$xyz_file" > "$gcrt_file"

antechamber -i "$gcrt_file" -fi gcrt \
    -o vg30-6-cl_gaff2.mol2 -fo mol2 \
    -c bcc -nc 1 -m 1 -rn VG3 -at gaff2

parmchk2 -i vg30-6-cl_gaff2.mol2 -f mol2 \
    -o vg30-6-cl_gaff2.frcmod -s gaff2

antechamber -i vg30-6-cl_gaff2.mol2 -fi mol2 \
    -o vg30-6-cl.pdb -fo pdb
