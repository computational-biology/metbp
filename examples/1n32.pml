load 1n32.cif
select polyall, polymer
select waterall, solvent
hide everything,1n32
set cartoon_ring_mode, 0, polyall
set cartoon_ladder_mode, 0, polyall
show cartoon, polyall
select MG1563A, (resi 1563 and chain A)
set sphere_scale, 0.40, MG1563A
show spheres, MG1563A
select tmp, (resi    858 and chain A  )  (resi    828 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 858 and chain A and name N7)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1563 and chain A),  (resi 858 and chain A and name N7), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp,  (resi 869 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color yellow, tmp
util.cbay tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 869 and chain A and name N7)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1563 and chain A),  (resi 869 and chain A and name N7), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1568A, (resi 1568 and chain A)
set sphere_scale, 0.40, MG1568A
show spheres, MG1568A
select tmp, (resi   1370 and chain A  )  (resi   1352 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1370 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1568 and chain A),  (resi 1370 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1576A, (resi 1576 and chain A)
set sphere_scale, 0.40, MG1576A
show spheres, MG1576A
select tmp, (resi    299 and chain A  )  (resi    566 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 299 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1576 and chain A),  (resi 299 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp,  (resi 558 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color green, tmp
util.cbag tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 558 and chain A and name OP1)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1576 and chain A),  (resi 558 and chain A and name OP1), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1577A, (resi 1577 and chain A)
set sphere_scale, 0.40, MG1577A
show spheres, MG1577A
select tmp, (resi    324 and chain A  )  (resi    109 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 324 and chain A and name N7)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1577 and chain A),  (resi 324 and chain A and name N7), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1586A, (resi 1586 and chain A)
set sphere_scale, 0.40, MG1586A
show spheres, MG1586A
select tmp, (resi   1511 and chain A  )  (resi   1524 and chain A  )  (resi    767 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1511 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1586 and chain A),  (resi 1511 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1595A, (resi 1595 and chain A)
set sphere_scale, 0.40, MG1595A
show spheres, MG1595A
select tmp, (resi    372 and chain A  )  (resi    389 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 372 and chain A and name O2)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1595 and chain A),  (resi 372 and chain A and name O2), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi    387 and chain A  )  (resi    376 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 387 and chain A and name O4)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1595 and chain A),  (resi 387 and chain A and name O4), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1598A, (resi 1598 and chain A)
set sphere_scale, 0.40, MG1598A
show spheres, MG1598A
select tmp, (resi   1199 and chain A  )  (resi   1058 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1199 and chain A and name O4)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1598 and chain A),  (resi 1199 and chain A and name O4), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1606A, (resi 1606 and chain A)
set sphere_scale, 0.40, MG1606A
show spheres, MG1606A
select tmp,  (resi 121 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color yellow, tmp
util.cbay tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 121 and chain A and name N3)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1606 and chain A),  (resi 121 and chain A and name N3), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi    124 and chain A  )  (resi    237 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 124 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1606 and chain A),  (resi 124 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi    236 and chain A  )  (resi    125 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 236 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1606 and chain A),  (resi 236 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1609A, (resi 1609 and chain A)
set sphere_scale, 0.40, MG1609A
show spheres, MG1609A
select tmp, (resi     98 and chain A  )  (resi     70 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 98 and chain A and name O4)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1609 and chain A),  (resi 98 and chain A and name O4), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1610A, (resi 1610 and chain A)
set sphere_scale, 0.40, MG1610A
show spheres, MG1610A
select tmp, (resi    105 and chain A  )  (resi     62 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 105 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1610 and chain A),  (resi 105 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1614A, (resi 1614 and chain A)
set sphere_scale, 0.40, MG1614A
show spheres, MG1614A
select tmp,  (resi 3 and chain Z)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color red, tmp
util.cbao tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 3 and chain Z and name O2')
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1614 and chain A),  (resi 3 and chain Z and name O2'), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1625A, (resi 1625 and chain A)
set sphere_scale, 0.40, MG1625A
show spheres, MG1625A
select tmp,  (resi 181 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color red, tmp
util.cbao tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 181 and chain A and name O2')
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1625 and chain A),  (resi 181 and chain A and name O2'), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp, (resi    183 and chain A  )  (resi    194 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 183 and chain A and name N7)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1625 and chain A),  (resi 183 and chain A and name N7), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1628A, (resi 1628 and chain A)
set sphere_scale, 0.40, MG1628A
show spheres, MG1628A
select tmp, (resi    758 and chain A  )  (resi    583 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 758 and chain A and name N7)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1628 and chain A),  (resi 758 and chain A and name N7), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1633A, (resi 1633 and chain A)
set sphere_scale, 0.40, MG1633A
show spheres, MG1633A
select tmp, (resi    516 and chain A  )  (resi    533 and chain A  )  (resi    519 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 516 and chain A and name O4)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1633 and chain A),  (resi 516 and chain A and name O4), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp,  (resi 533 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color green, tmp
util.cbag tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 533 and chain A and name OP1)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1633 and chain A),  (resi 533 and chain A and name OP1), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1642A, (resi 1642 and chain A)
set sphere_scale, 0.40, MG1642A
show spheres, MG1642A
select tmp, (resi   1486 and chain A  )  (resi   1414 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1486 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1642 and chain A),  (resi 1486 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1643A, (resi 1643 and chain A)
set sphere_scale, 0.40, MG1643A
show spheres, MG1643A
select tmp, (resi   1416 and chain A  )  (resi   1484 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1416 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1643 and chain A),  (resi 1416 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1646A, (resi 1646 and chain A)
set sphere_scale, 0.40, MG1646A
show spheres, MG1646A
select tmp, (resi    786 and chain A  )  (resi    796 and chain A  )  (resi    696 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 786 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1646 and chain A),  (resi 786 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG441A, (resi 441 and chain A)
set sphere_scale, 0.40, MG441A
show spheres, MG441A
select tmp, (resi    853 and chain A  )  (resi    833 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 853 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 441 and chain A),  (resi 853 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1652A, (resi 1652 and chain A)
set sphere_scale, 0.40, MG1652A
show spheres, MG1652A
select tmp,  (resi 115 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color red, tmp
util.cbao tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 115 and chain A and name O2')
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1652 and chain A),  (resi 115 and chain A and name O2'), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select tmp,  (resi 116 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color green, tmp
util.cbag tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 116 and chain A and name OP2)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1652 and chain A),  (resi 116 and chain A and name OP2), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG465A, (resi 465 and chain A)
set sphere_scale, 0.40, MG465A
show spheres, MG465A
select tmp, (resi    229 and chain A  )  (resi    133 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 229 and chain A and name O4)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 465 and chain A),  (resi 229 and chain A and name O4), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1672A, (resi 1672 and chain A)
set sphere_scale, 0.40, MG1672A
show spheres, MG1672A
select tmp,  (resi 1048 and chain A)
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
color red, tmp
util.cbao tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1048 and chain A and name O2')
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1672 and chain A),  (resi 1048 and chain A and name O2'), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1674A, (resi 1674 and chain A)
set sphere_scale, 0.40, MG1674A
show spheres, MG1674A
select tmp, (resi   1417 and chain A  )  (resi   1483 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 1417 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1674 and chain A),  (resi 1417 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
select MG1676A, (resi 1676 and chain A)
set sphere_scale, 0.40, MG1676A
show spheres, MG1676A
select tmp, (resi     46 and chain A  )  (resi    395 and chain A  ) 
color white, tmp
set cartoon_ring_mode, 2, tmp
set cartoon_ladder_mode, 0, tmp
util.cbaw tmp
show_as cartoon, tmp
show line, tmp
select ligand, (resi 46 and chain A and name O6)
set sphere_scale, 0.40, ligand
show spheres,  ligand
distance d1,  (resi 1676 and chain A),  (resi 46 and chain A and name O6), 15, 4
set dash_width, 1.0, d1
set dash_gap, 0.2, d1
hide label, d1
