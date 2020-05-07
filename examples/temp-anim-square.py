import gmsh
import sys

dt = 10

gmsh.initialize(sys.argv)
gmsh.merge("temp-square.msh")

n_steps = int(gmsh.option.getNumber("View[0].NbTimeStep"))
times = []
temps = []

for i in range(n_steps):
  kind, tags, temp, t, _ = gmsh.view.getModelData(0, i)
  temps.append(temp)
  times.append(t)

end_time = t
temp_inst = [0] * len(temp)
view_inst = gmsh.view.add("temp_temp") # as in temperature and temporal

t = 0
step = 0
i = 1
while t < end_time:
  if t > times[i]:
    while times[i] < t:  
      i += 1
  alpha = (t-times[i-1])/(times[i]-times[i-1])  
  print(t,i,alpha)
  
  for j in range(len(temps[i])):
    temp_inst[j] = temps[i-1][j] + alpha * (temps[i][j] - temps[i-1][j])

  gmsh.view.addModelData(view_inst, step, "", kind, tags, temp_inst, t)
  
  step += 1
  t += dt

gmsh.view.remove(0)
gmsh.fltk.initialize()

for i in range(step):
  print(i)  
  gmsh.option.setNumber("View[0].TimeStep", i)
  gmsh.fltk.update()
  gmsh.write("temp-square-smooth-%03d.png" % i)


gmsh.finalize()
print("all frames dumped, now run")
print("ffmpeg -y -framerate 20 -f image2 -i temp-square-smooth-%03d.png temp-square-smooth.mp4")
print("to get a video")


