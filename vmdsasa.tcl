set outfile [open sasa.txt w]
set nf [molinfo top get numframes]
#mol new 60FF_ff14ipq_rep1.nowater.prmtop
#mol addfile combinedtraj0-300.dcd first 0 step 1 waitfor all
set all [atomselect top "resid 1 to 120"]
for {set i 0} {$i < $nf} {incr i} {
	$all frame $i
	$all update
	set sasa [measure sasa 1.4 $all]
	puts $outfile "Frame $i, SASA $sasa" 
}
close $outfile
