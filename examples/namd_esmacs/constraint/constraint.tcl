mol load pdb /home/dave/tmp/bac-store/default/egfr-a698g/com/drug/aee/2j6m-a698g/build/complex.pdb
set a [atomselect top all]

set outfile [open tmp_cbv w]
set minmax [measure minmax $a]
set boxsize [vecsub [lindex $minmax 1] [lindex $minmax 0]]
set centre [measure center $a]

puts $outfile $boxsize
puts $outfile $centre

close $outfile

$a set beta 0
$a set occupancy 0
set b [atomselect top "(resid 1 to 288) and noh"]
$b set beta 4
set z [atomselect top "(not (resid 1 to 288) and not water and not type IM) and noh"]
$z set beta 4

set c [atomselect top "same residue as((resid 1 to 288) or (not (resid 1 to 288) and not water and not type IM)) and within 5 of resid 698 698"]
$c set beta 0

$a writepdb /home/dave/tmp/bac-store/default/egfr-a698g/com/drug/aee/2j6m-a698g/constraint/f4.pdb

quit
