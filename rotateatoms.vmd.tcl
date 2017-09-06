# Rotates a group of atoms (in the form of an atomselect), useful for
# positioning a ligand in its site
proc rotateatoms {m x y z} {
  set loc [measure center $m]
  $m moveby [vecinvert $loc]
  $m move [transaxis x $x]
  $m move [transaxis y $y]
  $m move [transaxis z $z]
  $m moveby $loc
}
