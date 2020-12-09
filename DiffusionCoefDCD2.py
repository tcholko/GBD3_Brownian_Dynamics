import sys
import numpy as np
import pytraj as pt

if len(sys.argv) < 6:
   sys.exit("Input error\nUsage: program.py traj.dcd traj.prmtop tstep(ps) tint(measure over this many frames) trest(frames between time origins)")
   
traj = pt.load(sys.argv[1], sys.argv[2])
tstep = float(sys.argv[3]); interval = int(sys.argv[4]); trest = int(sys.argv[5]);
com1 = 0; com2 = 0
tot_mass = 0
nframe = traj.n_frames
nrest = int((nframe-interval)/trest)
if nrest < 1:
   nrest = 1

print(str(nrest) + " restarts will be used")
natoms = traj.top.n_atoms
tdisp2 = 0

for m in range(natoms): tot_mass += traj.top[m].mass
print("Ligand mass is " + str(tot_mass) + " amu")

for i in range(nrest):
  for j in range(natoms):
    com1 += traj[i].atom(j) * traj.top[j].mass # Multiplies x y and z of atom j in frame i by mass all at once
    com2 += traj[i+interval].atom(j) * traj.top[j].mass 
  com1 /= tot_mass
  com2 /= tot_mass
  print("COM1 is " + str(com1) + " COM2 is " + str(com2))
  dx = com2[0]-com1[0]
  dy = com2[1]-com1[1]
  dz = com2[2]-com1[2]
  disp2 = (dx**2 + dy**2 + dz**2) 
  tdisp2 += disp2
 
avgdisp2 = tdisp2/nrest # Average over different time origins 
D = avgdisp2/(6*interval*tstep)
print(str(round(D*1e12*1e-16, 8)) + " cm2/s")
