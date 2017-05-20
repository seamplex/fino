rm -f *.pdf
rm -f *.vtk
rm -f *.dat
rm -f veeder.msh
rm -f cantilever-*.geo cantilever-*.msh
grep examples ../.gitignore | sed s/examples/\./ | xargs rm -f
