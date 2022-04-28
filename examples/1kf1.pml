load 1kf1.cif
select polyall, polymer
select waterall, solvent
hide everything,1kf1
set cartoon_ring_mode, 0, polyall
set cartoon_ladder_mode, 0, polyall
show cartoon, polyall
select K24A, (resi 24 and chain A)
set sphere_scale, 0.40, K24A
show spheres, K24A
select tmp, (resi      3 and chain A  )  (resi      9 and chain A  )  (resi     21 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 3 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 3 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      4 and chain A  )  (resi     10 and chain A  )  (resi     22 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 4 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 4 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      9 and chain A  )  (resi      3 and chain A  )  (resi     15 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 9 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 9 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     10 and chain A  )  (resi      4 and chain A  )  (resi     16 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 10 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 10 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     15 and chain A  )  (resi     21 and chain A  )  (resi      9 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 15 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 15 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     16 and chain A  )  (resi     22 and chain A  )  (resi     10 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 16 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 16 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     21 and chain A  )  (resi     15 and chain A  )  (resi      3 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 21 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 21 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     22 and chain A  )  (resi     16 and chain A  )  (resi      4 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 22 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 24 and chain A),  (resi 22 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select K25A, (resi 25 and chain A)
set sphere_scale, 0.40, K25A
show spheres, K25A
select tmp, (resi      2 and chain A  )  (resi      8 and chain A  )  (resi     20 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 2 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 2 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      3 and chain A  )  (resi      9 and chain A  )  (resi     21 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 3 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 3 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      8 and chain A  )  (resi      2 and chain A  )  (resi     14 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 8 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 8 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      9 and chain A  )  (resi      3 and chain A  )  (resi     15 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 9 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 9 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     14 and chain A  )  (resi     20 and chain A  )  (resi      8 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 14 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 14 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     15 and chain A  )  (resi     21 and chain A  )  (resi      9 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 15 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 15 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     20 and chain A  )  (resi     14 and chain A  )  (resi      2 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 20 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 20 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     21 and chain A  )  (resi     15 and chain A  )  (resi      3 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 21 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 25 and chain A),  (resi 21 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select K26A, (resi 26 and chain A)
set sphere_scale, 0.40, K26A
show spheres, K26A
select tmp, (resi      2 and chain A  )  (resi      8 and chain A  )  (resi     20 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 2 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 26 and chain A),  (resi 2 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi      8 and chain A  )  (resi      2 and chain A  )  (resi     14 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 8 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 26 and chain A),  (resi 8 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     14 and chain A  )  (resi     20 and chain A  )  (resi      8 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 14 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 26 and chain A),  (resi 14 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi     20 and chain A  )  (resi     14 and chain A  )  (resi      2 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 20 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 26 and chain A),  (resi 20 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
