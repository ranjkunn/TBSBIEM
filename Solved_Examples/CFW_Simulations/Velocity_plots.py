import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import os

# Parameters
mu_s = [0.3, 0.45, 0.5, 0.55, 0.677, 0.99]
v0_by_cs = [0.622, 0.715, 0.77, 0.802, 0.85, 0.865]

mu_s0 = 0.677
foldername = './'
resultfilename = 'mu677bg_'
magni_wave = 1
dia = 1000
Dia = 'D1000_'

# Determine v0 based on mu_s0
if mu_s0 == mu_s[0]:
    v0 = v0_by_cs[0]
elif mu_s0 == mu_s[1]:
    v0 = v0_by_cs[1]
elif mu_s0 == mu_s[2]:
    v0 = v0_by_cs[2]
elif mu_s0 == mu_s[3]:
    v0 = v0_by_cs[3]
elif mu_s0 == mu_s[4]:
    v0 = v0_by_cs[4]
elif mu_s0 == mu_s[5]:
    v0 = v0_by_cs[5]
else:
    raise ValueError('Try only 0.45 0.5 0.55 0.677 0.99')

n = 1626
Lx = 30000
Lz = 15000
dx = 15000/512  # 29.296875
dt = 5.5/1626 #3.3830109699769053E-003   # 0.003383011

# Load initial data
filename1 = os.path.join(foldername, 'dats/SlipFrontOut.dat')
t1 = np.loadtxt(filename1)

# First figure - Slip front plots
diff_incs = 296
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)

# Create asperity polygon
p = RegularPolygon((-11500/Lx+0.5, 0.5), 1000, radius=dia/Lx, fill=False)
ax1.add_patch(p)

h1_lines = []
h2_lines = []

j = 0
for i in range(296, n+1, diff_incs):
    try:
        filename2 = os.path.join(foldername, f'txts/SlipFrontOutNew{i:04d}.txt')
        t2 = np.loadtxt(filename2)
        j += 1
        
        # Plot cracktip and slip front
        h1 = ax1.axvline(t1[i, 1]/Lx + 0.5, color='k')
        h2, = ax1.plot(t2[1:-1, 0]/Lx + 0.5, t2[1:-1, 1]/Lz + 0.5, '-', markersize=1)
        
        h1_lines.append(h1)
        h2_lines.append(h2)
    except:
        continue

ax1.set_xlim([0, 1])
ax1.set_ylim([0, 1])
ax1.set_xlabel('Cracktip location (x/$L_x$)')
ax1.set_ylabel('z/$L_z$')

# Create legend
legend_labels = ['Cracktip no asperity', 'Asperity'] + \
               [f'Inc. = {i}' for i in range(296, n+1, diff_incs)[:len(h2_lines)]]
ax1.legend([h1_lines[0], p] + h2_lines, legend_labels, 
           loc='upper right', ncol=1, fontsize=10, framealpha=0)

plt.savefig(f'{resultfilename}{Dia}mus1pt1_mur0pt700_Slip.jpeg')
# plt.savefig(f'{resultfilename}{Dia}mus1pt1_mur0pt700_Slip.fig')

# Second figure - Velocity plots
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)

# Create asperity polygon again
p = RegularPolygon((-11500/Lx+0.5, 0.5), 1000, radius=dia/Lx, fill=False)
ax2.add_patch(p)

h1_vel_lines = []
h2_vel_lines = []

j = 0
for i in range(296, n+1, diff_incs):
    try:
        prev_i = i - diff_incs
        filename_prev = os.path.join(foldername, f'txts/SlipFrontOutNew{prev_i:04d}.txt')
        filename_curr = os.path.join(foldername, f'txts/SlipFrontOutNew{i:04d}.txt')
        
        t2 = np.loadtxt(filename_prev)
        t3 = np.loadtxt(filename_curr)
        
        t4 = np.empty((len(t3)-1, 2))
        t4[:, 0] = (t3[1:, 0] - t2[1:, 0])/(dt*diff_incs)
        t4[:, 1] = (t3[1:, 1] + t2[1:, 1])/2
        
        j += 1
        
        # Plot velocity data
        h1_vel = ax2.axvline(t1[i, 1]/Lx + 0.5, color='k')
        h2_vel, = ax2.plot((v0*3464*dt*i)/Lx + (t4[:, 0]/3464 - v0)/magni_wave, 
                          t4[:, 1]/Lz + 0.5, '-', markersize=1)
        
        h1_vel_lines.append(h1_vel)
        h2_vel_lines.append(h2_vel)
    except:
        continue

ax2.set_xlim([0, 1])
ax2.set_ylim([0, 1])
ax2.set_xlabel('(v-v_0)/c_s+ Const(v_0 t/L_z)')
ax2.set_ylabel('z/$L_z$')

# Create legend for velocity plot
ax2.legend([h1_vel_lines[0], p] + h2_vel_lines, legend_labels, 
           loc='upper right', ncol=1, fontsize=10, framealpha=0)

plt.savefig(f'{resultfilename}{Dia}mus1pt1_mur0pt700_Vel.jpeg')
# plt.savefig(f'{resultfilename}{Dia}mus1pt1_mur0pt700_Vel.fig')

# plt.show()