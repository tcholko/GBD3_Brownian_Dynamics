import math 
import sys

if len(sys.argv) < 5:
   sys.exit("Input error\nUsage: program.py traj.pdb tstep(ps) natoms nframes")

tstep = float(sys.argv[2]); natoms = int(sys.argv[3]); nframe = int(sys.argv[4]);
traj = sys.argv[1]
#nrest = int((nframe-tint)/trest)

xc = []; yc = []; zc = []
mass = []; tot_mass = 0
comx1 = 0; comy1 = 0; comz1 = 0
comx2 = 0; comy2 = 0; comz2 = 0
tdisp2 = 0
# Get coordinates and assign atomic masses
with open(traj, "r") as file:
       for line in file:
          if line.startswith("ATOM"):    
             xc.append(line[30:37])
             yc.append(line[38:45])
             zc.append(line[46:53])

          if line[13:14] == "H": mass.append(1.001)
          if line[13:14] == "C": mass.append(12.01)
          if line[13:14] == "O": mass.append(16.00)
          if line[13:14] == "N": mass.append(14.01)
          if line[13:14] == "P": mass.append(30.97)
          
for i in range (natoms): tot_mass += mass[i]
print("Ligand mass is " + str(tot_mass) + " amu")

for i in range(nframe-1):
    # Calculate ceneter-of-mass in frame i and i+1
    for j in range(natoms):
       comx1 += float(xc[j+i*natoms]) * float(mass[j+i*natoms]) 
       comy1 += float(yc[j+i*natoms]) * float(mass[j+i*natoms])
       comz1 += float(zc[j+i*natoms]) * float(mass[j+i*natoms]) 

       comx2 += float(xc[j+(i+1)*natoms]) * float(mass[j+(i+1)*natoms])
       comy2 += float(yc[j+(i+1)*natoms]) * float(mass[j+(i+1)*natoms])
       comz2 += float(zc[j+(i+1)*natoms]) * float(mass[j+(i+1)*natoms])

    comx1 /= tot_mass; comy1 /= tot_mass; comz1 /= tot_mass
    comx2 /= tot_mass; comy2 /= tot_mass; comz2 /= tot_mass
    #print("COM in frame " + str(i) + " = " + str(comx1) + " " + str(comy1) + " " + str(comz1))
    # Calculate displacement of center-of-mass
    x_disp = float(comx2 - comx1)
    y_disp = float(comy2 - comy1)
    z_disp = float(comz2 - comz1)
    disp2 = (x_disp**2 + y_disp**2 + z_disp**2)
    tdisp2 += disp2
    print(round(tdisp2, 3))
  
D = tdisp2/(2 * 3 * nframe*tstep) 
print("D = " + str(round(D*1e12*1e-16, 8)) + " cm^2/s")
