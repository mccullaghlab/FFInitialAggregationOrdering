set outfile [open sasa_backbone.txt w]
set nf [molinfo top get numframes]
#mol new 60FF_ff14ipq_rep1.nowater.prmtop
#mol addfile combinedtraj0-300.dcd first 0 step 1 waitfor all
gets stdin selmode
set sel [atomselect top "$selmode"]
#set not [atomselect top "resid 1 to 2 and sidechain"]
set all [atomselect top "protein"]
for {set i 0} {$i < $nf} {incr i} {
	$all frame $i
	$all update
	set sasa [measure sasa 1.4 $all -restrict $sel]
	puts $outfile "Frame $i, SASA $sasa" 
}
close $outfile
