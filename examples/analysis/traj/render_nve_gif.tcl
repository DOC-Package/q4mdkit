if {[llength $argv] < 3} {
    puts "Usage: vmd -dispdev text -e render_nve_gif.tcl -args structure.pdb trajectory.dcd output_dir ?frame_stride? ?zoom_scale? ?width? ?height?"
    exit 1
}

set structure_file [lindex $argv 0]
set trajectory_file [lindex $argv 1]
set output_dir [lindex $argv 2]
set frame_stride 4
set zoom_scale 1.8
set image_width 1600
set image_height 1600

if {[llength $argv] >= 4} {
    set frame_stride [lindex $argv 3]
}

if {[llength $argv] >= 5} {
    set zoom_scale [lindex $argv 4]
}

if {[llength $argv] >= 6} {
    set image_width [lindex $argv 5]
}

if {[llength $argv] >= 7} {
    set image_height [lindex $argv 6]
}

if {$frame_stride < 1} {
    puts "frame_stride must be >= 1"
    exit 1
}

if {$zoom_scale <= 0} {
    puts "zoom_scale must be > 0"
    exit 1
}

if {$image_width < 1 || $image_height < 1} {
    puts "image dimensions must be >= 1"
    exit 1
}

file mkdir $output_dir

mol new $structure_file type pdb waitfor all
mol addfile $trajectory_file type dcd waitfor all

set molid [molinfo top]
set numframes [molinfo $molid get numframes]

display resize $image_width $image_height
display projection Orthographic
axes location Off
color Display Background white

mol delrep 0 $molid
mol representation Licorice 0.2 18 18
mol color Name
mol selection all
mol material Opaque
mol addrep $molid

animate goto 0

set refsel [atomselect $molid all frame 0]
set refcenter [measure center $refsel weight mass]
$refsel moveby [vecscale -1.0 $refcenter]

display resetview
rotate x by -70
rotate y by -20
scale by $zoom_scale

for {set frame 0} {$frame < $numframes} {incr frame $frame_stride} {
    animate goto $frame
    set framesel [atomselect $molid all frame $frame]
    set framecenter [measure center $framesel weight mass]
    $framesel moveby [vecscale -1.0 $framecenter]
    set outfile [format "%s/frame_%04d.tga" $output_dir $frame]
    render TachyonInternal $outfile
    $framesel delete
}

$refsel delete

exit