package require psfgen
topology /home/dave/dev/BacBuilder/charmm_toppar/top_all27_prot_lipid.inp
alias residue ABA CYS
alias residue HIS HSE
alias residue HOH TIP3
alias atom HOH O OH2
alias atom ILE CD1 CD
alias residue NLE MET

segment A {
pdb chainA.pdb
mutate 698 GLY
}

coordpdb chainA.pdb A

guesscoord
writepdb tmp_mut.pdb
writepsf tmp_mut.psf
mol load psf tmp_mut.psf pdb tmp_mut.pdb
[atomselect top "segid A and noh"] writepdb chainA.pdb
exit
