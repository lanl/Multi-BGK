import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

base_filename = '../Data/cutout_allshoulder_yukawa_test34'

x = np.loadtxt(base_filename + '_x')
x_mm = x * 10;
Nx = x.size

x_start = int(Nx/4)
x_end = x_start + int(Nx/2)

#x_start = 0
#x_end = -1

Lv = 8648.171611
Nv = 60
Nv_slice = int(Nv/2)
v = np.linspace(-Lv, Lv, Nv)
v_ms = v / 100.


kb = 1.60218e-12 #erg / eV
erg_to_J = 1e-7 # J / erg
e0 = 1.602e-19 # C

#this is the dump number
step = 0

##this is how many cycles between dumps
t_offset = 50
dt = 1e-9


FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Distribution - Yuakwa', artist='Matplotlib',
                comment='Density plot')
writer = FFMpegWriter(fps=15, metadata=metadata, codec='mpeg2video')


fig = plt.figure()

num_ranks = 110

with writer.saving(fig, base_filename + ".mp4", 100):
    for k in np.arange(1,101):
        plt.cla()
        t_ms = t_offset*dt / 1e-6
        print("Time ", step*dt/1e-6)
        dist = np.fromfile(base_filename + '_spec0.step' + str(t_offset) + '_rank0.dat')
        nx_dist = int(dist.size / (Nv*Nv*Nv))
        dist_3d = np.reshape(dist, (nx_dist,Nv,Nv,Nv))
        
        for rank in range(1, num_ranks):
            dist_rank = np.fromfile(base_filename + '_spec0.step' + str(t_offset) + '_rank' + str(rank) + '.dat')
            nx_dist = int(dist_rank.size / (Nv*Nv*Nv))
            dist_rank_3d = np.reshape(dist_rank, (nx_dist,Nv,Nv,Nv))
            dist_3d = np.concatenate((dist_3d, dist_rank_3d))
            #print(dist_3d.shape)
            
        print(dist_3d.shape)
            
        #X, V = np.mgrid[x,v]
        plt.pcolor(x_mm[x_start:x_end], v_ms, np.transpose(dist_3d[x_start:x_end,:,Nv_slice,Nv_slice]))
        #plt.ylim([-Lv, Lv])
        plt.ylabel('Velocity [m/s]')
        plt.xlabel('Location [mm]')
        plt.title('Distribution at t = ' + str(t_ms) + ' microsec')
        plt.savefig(base_filename + str(t_ms) + 'usec_dist.png')
        writer.grab_frame()
