# VMD script to render a DCD trajectory and assemble it into a GIF.
# Usage:
#   vmd -e /path/to/make_gif_from_dcd.tcl
# Optional environment variables:
#   VMD_PDB_PATH       : topology PDB path
#   VMD_DCD_PATH       : trajectory DCD path
#   VMD_GIF_PATH       : output GIF path
#   VMD_FRAME_DIR      : directory for rendered frames
#   VMD_QM_SELECTION    : QM-region atom selection (default: beta > 50)
#   VMD_NEARBY_CUTOFF   : MM neighborhood cutoff in Angstrom (default: 8.0)
#   VMD_GIF_STRIDE     : render every Nth frame (default: 10)
#   VMD_GIF_MAX_FRAMES : maximum frames to render, 0 means all (default: 0)
#   VMD_GIF_DELAY      : ImageMagick delay in 1/100 s (default: 8)
#   VMD_KEEP_FRAMES    : 1 keeps rendered frames, 0 removes them after GIF creation
#   VMD_ALIGN_SEL      : atom selection used for alignment (default: beta > 50)
#   VMD_PBC_COMPOUND    : pbc wrap compound mode (default: fragment)
#   VMD_RENDERER       : snapshot or TachyonInternal (default: snapshot)

proc env_or_default {name default_value} {
    if {[info exists ::env($name)] && $::env($name) ne ""} {
        return $::env($name)
    }
    return $default_value
}

proc resolve_path {base_dir path_value} {
    if {[file pathtype $path_value] eq "absolute"} {
        return [file normalize $path_value]
    }
    return [file normalize [file join $base_dir $path_value]]
}

proc render_frame {renderer output_file} {
    if {$renderer eq "TachyonInternal"} {
        render TachyonInternal $output_file
    } else {
        render snapshot $output_file
    }
}

proc enable_dynamic_selection {molid repid} {
    if {[catch {mol selupdate $repid $molid on}]} {
        mol selupdate $repid $molid 1
    }
}

proc main {} {
    package require pbctools

    set script_dir [file dirname [file normalize [info script]]]
    set pdb_path [resolve_path $script_dir [env_or_default VMD_PDB_PATH "pentacene_qmmm.pdb"]]
    set dcd_path [resolve_path $script_dir [env_or_default VMD_DCD_PATH "../nve/output/nve.dcd"]]
    set gif_path [resolve_path $script_dir [env_or_default VMD_GIF_PATH "nve.gif"]]
    set frame_dir [resolve_path $script_dir [env_or_default VMD_FRAME_DIR "gif_frames"]]
    set qm_sel_text [env_or_default VMD_QM_SELECTION "beta > 50"]
    set nearby_cutoff [env_or_default VMD_NEARBY_CUTOFF "8.0"]
    set frame_stride [expr {int([env_or_default VMD_GIF_STRIDE "10"])}]
    set max_frames [expr {int([env_or_default VMD_GIF_MAX_FRAMES "0"])}]
    set gif_delay [expr {int([env_or_default VMD_GIF_DELAY "8"])}]
    set keep_frames [expr {int([env_or_default VMD_KEEP_FRAMES "0"])}]
    set align_sel_text [env_or_default VMD_ALIGN_SEL $qm_sel_text]
    set pbc_compound [env_or_default VMD_PBC_COMPOUND "fragment"]
    set renderer [env_or_default VMD_RENDERER "snapshot"]
    set nearby_mm_sel [format "same fragment as within %s of (%s) and not (%s)" $nearby_cutoff $qm_sel_text $qm_sel_text]

    if {![file exists $pdb_path]} {
        error "PDB file not found: $pdb_path"
    }
    if {![file exists $dcd_path]} {
        error "DCD file not found: $dcd_path"
    }
    if {$frame_stride < 1} {
        error "VMD_GIF_STRIDE must be >= 1"
    }

    file mkdir $frame_dir

    mol new $pdb_path type pdb waitfor all
    animate read dcd $dcd_path waitfor all top
    pbc wrap -center com -centersel $align_sel_text -compound $pbc_compound -all
    mol delrep 0 top

    mol representation Lines 1.000000
    mol color Beta
    mol selection $nearby_mm_sel
    mol material Opaque
    mol addrep top
    enable_dynamic_selection top 0

    mol representation CPK 1.000000 0.300000 12.000000 12.000000
    mol color Name
    mol selection $qm_sel_text
    mol material Opaque
    mol addrep top
    enable_dynamic_selection top 1

    display projection orthographic
    display depthcue off
    axes location off
    color Display Background white
    display resize 1000 1000

    set qm_sel [atomselect top $qm_sel_text frame 0]
    if {[$qm_sel num] > 0} {
        set center [measure center $qm_sel]
        molinfo top set center [list $center]
    }

    set num_frames [molinfo top get numframes]
    if {$num_frames < 1} {
        error "No trajectory frames were loaded from: $dcd_path"
    }

    if {[$qm_sel num] > 0} {
        set ref_sel [atomselect top $align_sel_text frame 0]
        set all_sel [atomselect top "all"]
        for {set frame 0} {$frame < $num_frames} {incr frame} {
            set move_sel [atomselect top $align_sel_text frame $frame]
            if {[$move_sel num] == [$ref_sel num]} {
                set transform [measure fit $move_sel $ref_sel]
                $all_sel frame $frame
                $all_sel move $transform
            }
            $move_sel delete
        }
        $all_sel delete
        $ref_sel delete
    }

    puts "QM selection: $qm_sel_text"
    puts "MM neighborhood selection: $nearby_mm_sel"

    set rendered_files {}
    set rendered_count 0
    for {set frame 0} {$frame < $num_frames} {incr frame $frame_stride} {
        if {$max_frames > 0 && $rendered_count >= $max_frames} {
            break
        }
        animate goto $frame
        set frame_file [file join $frame_dir [format "frame_%05d.tga" $rendered_count]]
        render_frame $renderer $frame_file
        lappend rendered_files $frame_file
        incr rendered_count
        puts [format "Rendered frame %d/%d -> %s" [expr {$frame + 1}] $num_frames $frame_file]
    }

    if {$rendered_count == 0} {
        error "No frames were rendered. Check VMD_GIF_STRIDE and VMD_GIF_MAX_FRAMES."
    }

    set convert_cmd [auto_execok convert]
    if {$convert_cmd eq ""} {
        puts "ImageMagick 'convert' was not found. Rendered frames were kept in: $frame_dir"
        puts "Create the GIF manually, for example: convert -delay $gif_delay -loop 0 $frame_dir/frame_*.tga $gif_path"
        return
    }

    eval exec [list $convert_cmd -delay $gif_delay -loop 0] $rendered_files [list $gif_path]
    puts "GIF created: $gif_path"

    if {!$keep_frames} {
        foreach frame_file $rendered_files {
            file delete -force $frame_file
        }
        if {[llength [glob -nocomplain [file join $frame_dir *]]] == 0} {
            file delete -force $frame_dir
        }
        puts "Temporary frames removed"
    } else {
        puts "Rendered frames kept in: $frame_dir"
    }
}

if {[catch {main} err]} {
    puts stderr $err
    exit 1
}

exit