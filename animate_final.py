import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

pre1 = 'rw_spiral25_in4True_sideview_drag_p2_t'
pre2 = '_w5_kim'
pre3 = '.npy'

template = np.load('rw_spiral25_in4True_sideview_drag_p2_t1_w5_kim0.npy')
shape = np.shape(template)
grid = np.zeros((710,shape[0],shape[1]),dtype='float')

for i in range(50):
  print i
  for j in range(710):
      x = np.load(pre1 + str(j) + pre2 + str(i) + pre3).astype('float')
      grid[j] += x

grid /= 50

grid2 = np.swapaxes(grid,1,2)
grid2 = grid2.tolist()
for i in grid2:
    for j in i:
        j += [0] * 300
grid2 = np.array(grid2)
grid2 = np.swapaxes(grid2,1,2)

profile = np.load('rw_spiral25_in4True_sideview_drag_p2_t709_w5_kimprofile.npy').astype(float)
profile /= 50.
profile = profile.astype('int').astype('float')

for i in profile:
    i -= 200.

profile = profile.astype('int')

for i in range(len(grid2)):
    for j in range(len(grid2[i][0])):
        grid2[i][:,j] = np.roll(grid2[i][:,j],profile[i][j])
grid3 = grid2[:,150:450,:]

fig = plt.figure(figsize = (27,5),dpi=120)
#
#
im = plt.imshow(grid3[0],cmap ='seismic',origin='lower',interpolation='none',aspect = 'auto')
plt.clim((0.,2.))
plt.title('Time = 0',fontsize = 24)
#
#
def updatefig(j):
#     # set the data in the axesimage object
    im.set_array(grid3[j])
    plt.title('Time = ' + str(j),fontsize=24)
#     # return the artists set
    return [im]
#
ani = animation.FuncAnimation(fig, updatefig, frames=range(710),
                              interval=50, blit=True)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=10000)
ani.save('im_finalp2_r50_single.mp4', writer=writer)
plt.show()
