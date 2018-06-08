import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

excess_net = []
for itt in range(50):
    print itt
    pre1 = 'rw2_spiral50_in4True_sideview_drag_p2_t'
    pre2 = '_w5_kim'
    pre3 = '.npy'
    # h = 400
    r = range(0,710)
    dump = []
    #
    # for i in r:
    #     x = np.load(pre1+str(i)+pre2+'0.npy')
    #     y = np.zeros((400-len(x[:,0]),300))
    #     z = np.vstack((x,y))
    #     dump.append(z)
    #     # print(np.shape(x))
    # dump = np.array(dump)
    #
    # fig = plt.figure(figsize = (28,8),dpi=120)
    small = []
    large = []

    for i in r:
        dump.append(np.load(pre1+str(i) +pre2 + str(itt) + pre3))

    for i in range(len(dump)):
        j = np.zeros((400,300))
        j[:len(dump[i])] = dump[i]
        dump[i] = j#[40:340]
    dump = np.array(dump)

    for i in dump:
        x = i[150:].flatten()
        small.append(len(x[x==1.]))
        large.append(len(x[x==2.]))

    small = np.array(small).astype('float')
    large = np.array(large).astype('float')
    large -= large[0]
    # large /= 2
    large_fraction = large / (large + small)
    excess = large_fraction# / (4./12.)
    excess_net.append(excess)

# plt.plot(excess)
# plt.show()


# im = plt.imshow(dump[0],cmap ='seismic',origin='lower',interpolation='none',aspect = 'auto')
# plt.clim((0.,2.))
# plt.title('Time = 0',fontsize = 24)
#
#
# def updatefig(j):
#     # set the data in the axesimage object
#     im.set_array(dump[j])
#     plt.title('Time = ' + str(j),fontsize=24)
#     # return the artists set
#     return [im]
#
# ani = animation.FuncAnimation(fig, updatefig, frames=range(700),
#                               interval=50, blit=True)
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=10000)
# ani.save('im_p2.mp4', writer=writer)
# plt.show()
