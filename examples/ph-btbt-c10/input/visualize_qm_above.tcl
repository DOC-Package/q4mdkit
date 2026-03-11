# VMD visualization script for QM molecule and adjacent layer above
# Molecules: QM(76), 30, 161, 2, 3 (layer z ≈ 75)
# Sorted by minimum atom-atom distance
# Usage: vmd -e visualize_qm_above.tcl

# Load structure
mol new /home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/ph-btbt-c10.pdb type pdb waitfor all

# Delete default representation
mol delrep 0 top

# Molecule 30 (min_dist=2.65Å): Licorice (red, ColorID 1)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 1
mol selection {index 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983}
mol material Opaque
mol addrep top

# Molecule 161 (min_dist=2.82Å): Licorice (orange, ColorID 3)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 3
mol selection {index 10304 10305 10306 10307 10308 10309 10310 10311 10312 10313 10314 10315 10316 10317 10318 10319 10320 10321 10322 10323 10324 10325 10326 10327 10328 10329 10330 10331 10332 10333 10334 10335 10336 10337 10338 10339 10340 10341 10342 10343 10344 10345 10346 10347 10348 10349 10350 10351 10352 10353 10354 10355 10356 10357 10358 10359 10360 10361 10362 10363 10364 10365 10366 10367}
mol material Opaque
mol addrep top

# Molecule 2 (min_dist=3.16Å): Licorice (green, ColorID 7)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 7
mol selection {index 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191}
mol material Opaque
mol addrep top

# Molecule 3 (min_dist=3.16Å): Licorice (cyan, ColorID 10)
mol representation Licorice 0.2 12.0 12.0
mol color ColorID 10
mol selection {index 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255}
mol material Opaque
mol addrep top

# QM molecule 76: CPK representation (element colors)
mol representation CPK 1.0 0.3 12.0 12.0
mol color Name
mol selection {index 4864 4865 4866 4867 4868 4869 4870 4871 4872 4873 4874 4875 4876 4877 4878 4879 4880 4881 4882 4883 4884 4885 4886 4887 4888 4889 4890 4891 4892 4893 4894 4895 4896 4897 4898 4899 4900 4901 4902 4903 4904 4905 4906 4907 4908 4909 4910 4911 4912 4913 4914 4915 4916 4917 4918 4919 4920 4921 4922 4923 4924 4925 4926 4927}
mol material Opaque
mol addrep top

# Change carbon color to gray for QM region
color Name C gray

# Display settings
display projection orthographic
display depthcue off
axes location off
color Display Background white

# Center on all selected molecules
set sel [atomselect top "index 128 to 255 1920 to 1983 4864 to 4927 10304 to 10367"]
set center [measure center $sel]
molinfo top set center [list $center]

# ============================================================
# Viewpoint save/restore functions
# ============================================================
proc save_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_above.tcl"}} {
    set fp [open $filename w]
    puts $fp "# VMD viewpoint settings - auto-generated"
    puts $fp "molinfo top set {rotate_matrix} {[molinfo top get rotate_matrix]}"
    puts $fp "molinfo top set {center_matrix} {[molinfo top get center_matrix]}"
    puts $fp "molinfo top set {scale_matrix} {[molinfo top get scale_matrix]}"
    puts $fp "molinfo top set {global_matrix} {[molinfo top get global_matrix]}"
    close $fp
    puts "Viewpoint saved to: $filename"
}

proc load_viewpoint {{filename "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_above.tcl"}} {
    if {[file exists $filename]} {
        source $filename
        puts "Viewpoint loaded from: $filename"
    } else {
        puts "Viewpoint file not found: $filename"
    }
}

# Auto-load viewpoint if exists
if {[file exists "/home/takahashi/python/q4mdkit/examples/ph-btbt-c10/input/viewpoint_qm_above.tcl"]} {
    load_viewpoint
}

puts ""
puts "=== QM + Adjacent Layer Above (z ≈ 75) ==="
puts "  Sorted by minimum atom-atom distance"
puts "  30: red (min_dist=2.65Å)"
puts "  161: orange (min_dist=2.82Å)"
puts "  2: green (min_dist=3.16Å)"
puts "  3: cyan (min_dist=3.16Å)"
puts "  76: gray (QM molecule)"
puts ""
puts "Toggle representations:"
puts "  mol showrep top 0 off/on  - Molecule 30 (red)"
puts "  mol showrep top 1 off/on  - Molecule 161 (orange)"
puts "  mol showrep top 2 off/on  - Molecule 2 (green)"
puts "  mol showrep top 3 off/on  - Molecule 3 (cyan)"
puts "  mol showrep top 4 off/on  - Molecule 76 (QM)"
puts ""
puts "Commands:"
puts "  save_viewpoint - Save current view"
puts "  load_viewpoint - Load saved view"
