rm -f *.png *.pdf
rm -f *.vtk
rm -f *.dat
rm -f veeder.msh
grep examples ../.gitignore | sed s/examples/\./ | xargs rm -f
