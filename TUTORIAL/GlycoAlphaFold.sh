#!/bin/bash
#(C) Mateusz Sikora 2021
### Get AlphaFold entry and run GlycoSHIELD on it, produces a trajectory with the largest number of conformers present on each glycosite.
### Example:
### ./GlycoAlphaFold.sh P13591
### where P13591 is a uniptor ID.

alphafold=$1
uniprot=${alphafold##*/}
pdbfile="https://alphafold.ebi.ac.uk/files/AF-${uniprot}-F1-model_v1.pdb"
infile=AF-${uniprot}-F1-model_v1.pdb

wget $pdbfile -O $infile 2>/dev/null

if [ -s "$infile" ]
then 
   :
else
   echo "AlphaFold entry does not exist or pdb file not found! Exiting..."
   exit 1
fi


cat <<EOT>aaa1
A XXX,YYY,ZZZ 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc A_YYY.pdb A_YYY.xtc
EOT

wget https://www.uniprot.org/uniprot/${uniprot} -q -O - | awk '{if($0~/Amino acid modifications/){ok++;next}; if(ok){print}}' | grep Glyco | python -c '
from __future__ import print_function
import sys
k=""
for line in sys.stdin:
   k+=line.strip()
for glycosite in k.split("key=Glycosylation\">")[1:]:
   print(glycosite.split("</a")[0])
' | while read line;  do xxx=$((${line}-1)); zzz=$((${line}+1)) ; sed 's/YYY/'${line}'/g' aaa1 | sed 's/XXX/'${xxx}'/g'|sed 's/ZZZ/'${zzz}'/g' ; done > input

if [ -s "$input" ]
then 
   :
else
   echo "UNIPROT entry does not contain any glycosites, exiting..."
   exit 1
fi

out=$(python GlycoSHIELD-0.1.py --inputfile input --protpdb $infile --mode CG --threshold 3.5)
smallestcommon=$(echo $out | tr ' ' '\n' | sort -n | awk '{if(NR==2){print $1}}')
python GlycoTRAJ-0.1.py --maxframe $smallestcommon --pdblist $(ls -v A_*.pdb | tr '\n' '\,' | sed 's/,*$//g') --xtclist $(ls -v A_*.xtc | tr '\n' '\,' | sed 's/,*$//g') --chainlist $(cat input | while read line; do echo "A"; done | tr '\n' '\,' | sed 's/,*$//g') --reslist $(cat input | while read line; do echo $line | awk '{split($2,a,","); print a[2]}'; done | tr '\n' '\,' | sed 's/,*$//g') --outname traj_$smallestcommon
