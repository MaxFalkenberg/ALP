import numpy as np
import matplotlib.pyplot as plt
import copy
from random import shuffle

def cd(v,a):
    vis = 1.51e-5
    Re = float(2. * a * v)/vis
    cd_raw = ((24./Re) * (1. + (0.27 * Re)) **  0.43) + (0.47 * (1. - np.exp(-0.04 * Re ** 0.38)))
    # cd = 0.5 + (2.* cd_raw)
    return cd_raw

def tra(v0,theta,a,dt,n):
    rho = 2650.
    m = (4./3.) * np.pi * (a ** 3) * rho
    x = np.zeros(n,dtype='float')
    z = np.zeros(n,dtype='float')
    v = np.zeros(n,dtype='float')
    vx = np.zeros(n,dtype='float')
    vz = np.zeros(n,dtype='float')
    v[0] = v0
    vx[0] = v0 * np.cos(theta)
    vz[0] = v0 * np.sin(theta)
    for i in range(1,n):
        dvx = drag_force(v[i-1],a) * np.absolute(np.cos(theta)) * dt * np.sign(vx[i-1]) / m
        dvz = drag_force(v[i-1],a) * np.absolute(np.sin(theta)) * dt * np.sign(vz[i-1]) / m
        # print vx[i-1], dvx, theta
        # print vz[i-1],dvz
        vx[i] = vx[i-1] - dvx
        vz[i] = vz[i-1] - (9.81 * dt) - dvz
        print (9.81 * dt) + dvz
        v[i] = ((float(vx[i]) ** 2.) + (float(vz[i]) ** 2)) ** 0.5
        theta = np.arctan(vz[i]/vx[i])
        x[i] = x[i-1] + vx[i] * dt
        z[i] = z[i-1] + vz[i] * dt
    return x,z


def drag_force(v,a):
    c = cd(v,a)
    rho_air = 1.2041
    f = 0.5 * rho_air * np.pi * (a**2) * (v**2) * c
    return f

def dv(vx,vz,dt,a):
    v = ((vx**2.)+(vz**2.))**0.5
    # print(vx,vz,v)
    c = cd(v,a)
    # print(c)
    coeff = (3. * c)/(8. * a)
    dvz = (-9.81 - ((coeff * (v**2)))) * dt * np.cos(np.arctan(vz/vx))
    dvx = (- (coeff * (v**2)))*dt * np.sin(np.arctan(vz/vx))
    return dvx,dvz

def traj_exact(v0,theta,a,dt,n):
    x = np.zeros(n)
    z = np.zeros(n)
    vx0 = np.cos(theta) * v0
    vz0 = np.sin(theta) * v0
    # print(vx0,vz0)
    vx = np.zeros(n)
    vz = np.zeros(n)
    vx[0] = vx0
    vz[0] = vz0
    # print(vx[:10],vz[:10])
    for i in range(1,n):
        dvx,dvz = dv(vx[i-1],vz[i-1],dt,a)
        vx[i] = vx[i-1] + dvx
        vz[i] = vz[i-1] + dvz
        x[i] = x[i-1] + vx[i-1] * dt
        z[i] = z[i-1] + vz[i-1] * dt
    return x,z

def results():
    # prange = [0.1,0.3,0.5,0.7,0.9]
    prange = [0.1,0.3,0.5]
    for n in range(10):
        print n
        for j in prange:
            a = Oslo(L = 250, w = 5, grain_throw = 'drag',p_large = j)
            a.run(1000)
            topview = np.vstack(a.topview)
            width = np.array(a.width_dump)
            depth = a.L - np.array(a.depth_dump)
            lr = np.array(a.lr)
            rwindow = np.array(a.removal_window_dump)
            a.removed_raw[0] = [a.p]
            removed_grains = np.array([np.mean(i) for i in a.removed_raw])-1
            p = str(a.p)
            p = p[2:]
            pre = 'kim' + 'L' + str(a.L) + 'w' + str(a.w) + 'pl' + p + a.grain_throw + 't' + str(a.t) + '_'
            np.save(pre + 'top'+str(n),topview)
            np.save(pre + 'width'+str(n),width)
            np.save(pre + 'depth'+str(n),depth)
            np.save(pre + 'lr'+str(n),lr)
            np.save(pre + 'rwindow'+str(n),rwindow)
            np.save(pre + 'rgrain'+str(n),removed_grains)

def sideview(rw = 5,largeprob = 0.2,spiral_on = True, throwmode = 'drag', gridL = 300,itts = 50):
    p0 =  largeprob
    mode = throwmode
    w = rw
    saverange = range(0,710)
    print'Test'
    # profile = np.zeros((710,300)).astype('float')
    for i in range(1,itts):
        print(i)
        a = Oslo(L = gridL, w = w, grain_throw = mode,p_large = p0,spiral = spiral_on)
        pre = 'rw2_spiral_in4'+str(spiral_on)+'_sideview_'+mode+'_p'+str(p0)[-1]+'_t709_w'+str(w)+'_kim'
        k = 0
        for j in saverange:
            a.run(j-k)
            k = j
            x = a.plot_pile()
        #     prof = []
        #     for q in a.grid:
        #         prof.append(len(q))
        #     prof = np.array(prof)
        #     profile[j]+=prof
            pre = 'rw2_spiral50_in4'+str(spiral_on)+'_sideview_'+mode+'_p'+str(p0)[-1]+'_t'+str(a.t)+'_w'+str(w)+'_kim'
            np.save(pre+str(i) ,x[150:])
        # a.run(710)
        # x = a.plot_pile()
        # np.save('spiral'+str(spiral_on)+'_sideview_'+mode+'_p'+str(p0)[-1]+'_t'+str(a.t)+'_w'+str(w)+'_kim'+str(i) ,x)
        # a.run(200)
        # x = a.plot_pile()
        # np.save('spiral'+str(spiral_on)+'_sideview_'+mode+'_p'+str(p0)[-1]+'_t'+str(a.t)+'_w'+str(w)+'_kim'+str(i) ,x)

        # topview = np.vstack(a.topview)
        # width = np.array(a.width_dump)
        # depth = a.L - np.array(a.depth_dump)
        # lr = np.array(a.lr)
        # rwindow = np.array(a.removal_window_dump)
        # avalanches = np.array(a.s)
        # l = []
        # s = []
        # for j in rwindow:
        #     l.append(len(j[j==1]))
        #     s.append(len(j[j==0]))
        # l = np.array(l).astype('float')
        # s = np.array(s).astype('float')
        # l*=2
        # ratio = l/(l+s)
        # a.removed_raw[0] = [a.p]
        # removed_grains = np.array([np.mean(j) for j in a.removed_raw])-1
        # p = str(a.p)
        # p = p[2:]
        # # pre = 'kim' + 'L' + str(a.L) + 'w' + str(a.w) + 'pl' + p + a.grain_throw + 't' + str(a.t) + '_'
        # # pre = 'spiral'+str(spiral_on)+'_' +mode+'_p'+str(p0)[-1]+'_t'+str(a.t)+'_w'+str(w)+'_kim'
        # np.save(pre + 'top'+str(i),topview)
        # np.save(pre + 'width'+str(i),width)
        # np.save(pre + 'depth'+str(i),depth)
        # np.save(pre + 'lr'+str(i),lr)
        # np.save(pre + 'ratio'+str(i),ratio)
        # np.save(pre + 'avalanches'+str(i),avalanches)
    # np.save(pre + 'profile',profile)
        # np.save(pre + 'rgrain'+str(i),removed_grains)
    # plt.figure()
    # plt.plot(ratio)
    # plt.figure()
    # plt.plot(avalanches)
    # plt.figure()
    # plt.imshow(x[160:],origin='lower',interpolation='none')
    # plt.figure()
    # plt.plot(profile[-1])
    # plt.show()


def density(L = 150.):
    x = np.linspace(0.1,L,1000).astype('float')
    d2 = (x**2)/(2*(L - x))
    d3 = (x ** 3)/(12*(((L/2)**2) - ((x/2)**2)))
    plt.plot(x,d2,label = '2D Density')
    plt.plot(x,d3,label = '3D Density')
    plt.plot([50,50],[0.1,100],label = 'Width = 50',ls = '--')
    plt.xlabel('Pit Width [Grain Units]', fontsize = 18)
    plt.ylabel('Grain Density', fontsize = 18)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc = 'best',fontsize=18)
    plt.savefig('graindensity.png')
    plt.show()

def lining_count(p = 0.5):
    N = 25
    x = [2,4,6,8,10,15,20,30,50,100,200,300,500,1000]
    dx = [2,2,2,2,2,5,5,10,20,50,100,100,200,500]
    yraw = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for i in range(N):
        print i
        a = Oslo(351,p_large=p,mode='dig',w=5)
        for j in range(len(dx)):
            a.run(dx[j])
            a.count()
            yraw[j].append(np.mean(a.lining))
    y = np.array([np.mean(i) for i in yraw])
    yerr = np.array([np.std(i,ddof=1) for i in yraw])
    return x,y,yerr



class Oslo:

    def __init__(self,L,p_large=0.5,mode='dig',centre=None,w=1,grain_throw = 'bucket',spiral = False):
        self.L = int(L)
        self.z = np.zeros((2,L),dtype='int8')
        self.p = p_large
        self.rdrop = []
        self.ldrop = []
        self.removed = []
        self.w = w
        self.width_dump = []
        self.depth_dump = []
        self.topview = []
        self.removed_raw = []
        self.lr = []
        self.s = []
        self.t = 0
        self.removal_window = np.zeros((self.w,self.w),dtype = 'float')
        self.removal_window_dump = []
        self.grain_throw = grain_throw
        self.spiral = spiral
        if mode == 'dig':
            # self.grains = [[0]*L for i in range(L)]
            # self.grains = np.random.rand(L,L) + 1 + p_large
            # self.grains = self.grains.astype('int').tolist()
            self.grains = [[] for i in range(self.L)]
            self.l = []
            for i in range(self.L):
                l = 0
                while np.sum(self.grains[i]) < self.L:
                    x = int(np.random.rand() + p_large + 1)
                    l += x
                    self.grains[i].append(x)
                self.l.append(l)
            self.l = np.array(self.l)
            self.z[1][:-1] = self.l[:-1] - self.l[1:]
            self.z[0][1:] = self.l[1:] - self.l[:-1]
            self.z[1][-1] = self.l[-1] - self.L
            self.z[0][0] = self.l[0] - self.L
            self.zc = np.zeros(L,dtype='int')
            for i in range(L):
                self.zc[i] = (2*self.grains[i][-2]) - (1*self.grains[i][-1]) + np.random.randint(2,4)
        elif mode == 'build':
            self.grains = [[] for i in range(L)]
        else:
            raise ValueError('mode variable must be set to \'build\' or \'dig\'.')
        self.mode = mode
        if centre==None:
            self.centre = int(L - 1)/2
            if self.spiral:
                self.wlist = []
                self.w0 = self.centre
                # for i in range(25,0,-1):
                #     self.wlist.append(self.centre + i)
                #     self.wlist.append(self.centre - i)
                #     self.wlist.append(self.centre + i)
                #     self.wlist.append(self.centre - i)
                #     self.wlist.append(self.centre + i)
                #     self.wlist.append(self.centre - i)
                #     self.wlist.append(self.centre + i)
                #     self.wlist.append(self.centre - i)
                for i in range(50,0,-1):
                    self.wlist.append(self.centre - i)
                    self.wlist.append(self.centre + i)
                    self.wlist.append(self.centre - i)
                    self.wlist.append(self.centre + i)
                    self.wlist.append(self.centre - i)
                    self.wlist.append(self.centre + i)
                    self.wlist.append(self.centre - i)
                    self.wlist.append(self.centre + i)



        else:
            self.centre = centre
        self.check = np.zeros(self.L,dtype='bool')

    def plot_pile(self):
        grid = copy.deepcopy(self.grains)
        for i in range(self.L):
            t = []
            for j in grid[i]:
                if j == 1:
                    t += [1]
                else:
                    t += [2,2]
            grid[i] = t
        self.grid = copy.deepcopy(grid)
        l = 0
        for i in grid:
            if len(i) > l:
                l = len(i)
        for i in range(len(grid)):
            if len(grid[i]) != l:
                # grid[i] = [0] * (l - len(grid[i])) + grid[i]
                grid[i] += [0] * ((l - len(grid[i])))
        # plt.clf()
        # plt.figure()
        # plt.imshow(np.array(grid).T, origin='lower', interpolation = 'none', aspect = 'auto')
        # plt.colorbar()
        # plt.show()
        return np.array(grid).T
    # @profile
    def count(self):
        l = True
        r = True
        len_t = np.array([np.sum(i) for i in self.grains])
        rlim = self.L-1
        llim = 0
        lining = [self.grains[self.centre][-1]]
        centre = self.centre
        if self.spiral:
            if centre < self.w0:
                centre = self.w0
        for i in range(centre+1,self.L):
            if r:
                # if len_t[i] > len_t[i-1]:
                #      lining += self.grains[i][len_t[i-1]:]
                # else:
                #     lining.append(self.grains[i][-1])
                lining.append(self.grains[i][-1])
                if len_t[i] >= self.L:
                    if len_t[i] >= len_t[i+1]:
                        rlim = i
                        r = False
                # else:
                #     r = False
        if self.spiral:
            if centre > self.w0:
                centre = self.w0
        for i in range(centre-1,-1,-1):
            if l:
                # if len_t[i] > len_t[i+1]:
                #     lining += self.grains[i][len_t[i+1]:]
                # else:
                #     lining.append(self.grains[i][-1])
                lining.append(self.grains[i][-1])
                if len_t[i] >= self.L:
                    if len_t[i] >= len_t[i-1]:
                        llim = i
                        l = False
                # else:
                #     l = False

        w = rlim - llim
        self.width = w
        self.width_dump.append(w)
        self.lr.append([llim,rlim])
        self.depth_dump.append(len_t[self.centre])
        self.topview.append([j[-1] for j in self.grains])
        self.lining = lining
        self.removed_raw.append([])
        for i in self.grains[:llim+1]:
            self.removed_raw[-1] += i[self.L:]
        for i in self.grains[rlim:]:
            self.removed_raw[-1] += i[self.L:]
        # print w, llim, rlim, lining
    # @profile
    def run(self,N=1):
        if self.mode=='build':
            self.check_sites = [self.centre]
            self.check[self.centre] = True
            for i in range(N):
                self.t += 1
                if np.random.rand() > self.p:
                    self.grains[self.centre].append(1)
                else:
                    self.grains[self.centre].append(2)
                self.z[:,self.centre] += 1
                if self.centre != self.L - 1:
                    self.z[:,self.centre+1][0] -= 1
                if self.centre != 0:
                    self.z[:,self.centre-1][1] -=1
                self.micro_run()
        else:
            # if N>=self.L:
            #     raise ValueError('N cannot exceed L for mode=\'dig\'')
            self.check[:] = True
            self.check_sites = range(self.L)
            for i in range(N):
                if self.spiral:
                    if self.t < len(self.wlist):
                        # print(1)
                        self.centre = self.wlist[self.t]
                    else:
                        self.centre = self.w0
                        # print(2)
                self.t+=1
                if np.sum(self.grains[self.centre]) <= self.w:
                    print('Digging ended because lower limit has been reached.')
                    print('Digging ended after ' + str(self.t) + ' iterations.')
                    break
                self.count()
                self.stemp = 0
                self.remove()
                self.micro_run()
                self.s.append(self.stemp)
    # @profile
    def remove(self):

        r = np.arange(self.centre,self.centre + self.w)
        r -= int(self.w - 1)/2
        rem_temp = np.zeros((self.w,self.w),dtype = 'float')
        r_depth = np.zeros_like(r)
        if np.min(r) < 0 or np.max(r) >= self.L:
            raise ValueError('Removal window size and location overlaps with system boundaries.')
        else:
            temp_store_g = []
            temp_store_z = []
            for i in range(self.w):
                for j in r:
                    if r_depth[j - r[0]] < self.w - 1:
                        g = self.grains[j].pop()
                        self.z[0][j] -= g
                        self.z[1][j] -= g
                        self.z[0][j+1] += g
                        self.z[1][j-1] += g
                        rem_temp[i][j-r[0]] += (g - 1)
                        z = np.sum(self.grains[j])
                        temp_store_g.append(g)
                        temp_store_z.append(z)
                        r_depth[j-r[0]] += g
                    elif r_depth[j - r[0]] == self.w - 1:
                        if self.grains[j][-1] == 1 or np.random.rand() > 0.5:
                            g = self.grains[j].pop()
                            self.z[0][j] -= g
                            self.z[1][j] -= g
                            self.z[0][j+1] += g
                            self.z[1][j-1] += g
                            rem_temp[i][j-r[0]] += (g - 1)
                            z = np.sum(self.grains[j])
                            temp_store_g.append(g)
                            temp_store_z.append(z)
                            r_depth[j-r[0]] += g
                        else:
                            g = 0
                            rem_temp[i][j-r[0]] += (g - 1)
                            z = np.sum(self.grains[j])
                            temp_store_g.append(g)
                            temp_store_z.append(z)
                            r_depth[j-r[0]] = self.w
                    else:
                        g = 0
                        rem_temp[i][j-r[0]] += (g - 1)
                        z = np.sum(self.grains[j])
                        temp_store_g.append(g)
                        temp_store_z.append(z)

                    # if self.grain_throw != 'bucket':
                    #     p = self.redistribute(g,i,j,z)
            k = 0
            # x = self.plot_pile()
            # plt.imshow(x,origin='lower',interpolation='none')
            # plt.show()
            self.micro_run()
            # self.plot_pile()
            if self.grain_throw != 'bucket':
                for i in range(self.w):
                    for j in r:
                        p = self.redistribute(temp_store_g[k],i,j,temp_store_z[k] + 5)
                        k += 1
            # x = self.plot_pile()
            # plt.imshow(x,origin='lower',interpolation='none')
            # plt.show()
            # self.z[:,r[0]][0] -= self.w
            # self.z[:,r[-1]][1] -= self.w
            # if r[0] != 0:
            #     self.z[:,r[0]-1][1] += self.w
            # if r[-1] != self.L - 1:
            #     self.z[:,r[-1]+1][0] += self.w
            self.removal_window += rem_temp
            self.removal_window_dump.append(rem_temp)
    # @profile
    def micro_run(self):
        s_total = 0
        cont = 1
        x = 0
        while cont:
            x += 1
            # if x >= 1000:
            #     self.plot_pile()
            sm = 0
            shuffle(self.check_sites)
            for i in self.check_sites:
                if len(self.grains[i]) < 2:
                    pass
                else:
                    zc = (2*self.grains[i][-2]) - (1*self.grains[i][-1]) + np.random.randint(2,4)
                    # if zc == 0:
                    #     zc = 2
                    # else:
                    #     zc += np.random.randint(2,4)
                    # if i == 0:
                    #     pass
                    # elif i == self.L -1:
                    #     pass
                    # else:
                    # lr = self.z[:,i]>zc
                    lr = self.z[:,i]>self.zc[i]
                    s = np.sum(lr)
                    sm += s
                    if s == 0:
                        pass
                    else:
                        s_total += 1
                        if s == 2:
                            if np.random.rand() > 0.5:
                                lr[0] = True
                                lr[1] = False
                            else:
                                lr[1] = True
                                lr[0] = False
                        if lr[0]:
                            if i == 0:
                                gsize = self.grains[i].pop()
                                self.ldrop.append(gsize)
                                self.z[:,i] -= 1 * gsize
                                self.z[:,i+1][0] += 1 *gsize
                            elif i == self.L - 1:
                                gsize = self.grains[i].pop()
                                self.grains[i-1].append(gsize)
                                self.z[:,i][0] -= 2 * gsize
                                self.z[:,i][1] -= 1 * gsize
                                self.z[:,i-1][0] += 1 * gsize
                                self.z[:,i-1][1] += 2 * gsize
                                self.z[:,i-2][1] -= 1 * gsize
                                if self.check[i-1] == False:
                                    self.check_sites.append(i-1)
                                    self.check[i-1] = True
                                self.zc[i-1] = (2*self.grains[i-1][-2]) - (1*self.grains[i-1][-1]) + np.random.randint(2,4)
                            else:
                                gsize = self.grains[i].pop()
                                self.grains[i-1].append(gsize)
                                self.z[:,i][0] -= 2 * gsize
                                self.z[:,i][1] -= 1 * gsize
                                self.z[:,i-1][0] += 1 * gsize
                                self.z[:,i-1][1] += 2 * gsize
                                self.z[:,i+1][0] += 1 * gsize
                                if i != 1:
                                    self.z[:,i-2][1] -=1 * gsize
                                if self.check[i-1] == False:
                                    self.check_sites.append(i-1)
                                    self.check[i-1] = True
                                self.zc[i-1] = (2*self.grains[i-1][-2]) - (1*self.grains[i-1][-1]) + np.random.randint(2,4)
                        elif lr[1]:
                            if i == 0:
                                gsize = self.grains[i].pop()
                                self.grains[i+1].append(gsize)
                                self.z[:,i][0] -= 1 * gsize
                                self.z[:,i][1] -= 2 * gsize
                                self.z[:,i+1][0] += 2 * gsize
                                self.z[:,i+1][1] += 1 * gsize
                                self.z[:,i+2][0] -= 1 * gsize
                                if self.check[i+1] == False:
                                    self.check_sites.append(i+1)
                                    self.check[i+1] = True
                                self.zc[i+1] = (2*self.grains[i+1][-2]) - (1*self.grains[i+1][-1]) + np.random.randint(2,4)
                            elif i == self.L - 1:
                                gsize = self.grains[i].pop()
                                self.rdrop.append(gsize)
                                self.z[:,i] -= 1 * gsize
                                self.z[:,i-1][1] += 1 * gsize
                            else:
                                gsize = self.grains[i].pop()
                                self.grains[i+1].append(gsize)
                                self.z[:,i][0] -= 1 * gsize
                                self.z[:,i][1] -= 2 * gsize
                                self.z[:,i+1][0] += 2 * gsize
                                self.z[:,i+1][1] += 1 * gsize
                                self.z[:,i-1][1] += 1 * gsize
                                if i != self.L - 2:
                                    self.z[:,i+2][0] -=1 * gsize
                                if self.check[i+1] == False:
                                    self.check_sites.append(i+1)
                                    self.check[i+1] = True
                                self.zc[i+1] = (2*self.grains[i+1][-2]) - (1*self.grains[i+1][-1]) + np.random.randint(2,4)
                        self.zc[i] = (2*self.grains[i][-2]) - (1*self.grains[i][-1]) + np.random.randint(2,4)
            if sm == 0:
                cont = 0
        self.stemp += s_total
    # @profile
    def redistribute(self,gsize,layer,x0,z0):
        if gsize == 0:
            pass
        else:
            # v0 = [70.,70.,70.,70.,70.]
            v0 = 70.
            dt = (np.random.rand()*20)-10
            dv = (np.random.rand()*60)-30
            # x,z = self.traj(gsize, v0[layer] + dv, 50. + dt)
            x,z = self.traj(gsize, v0 + dv, 50. + dt)
            # x,z = self.tra(gsize, (v0 + dv)/100, 50. + dt)
            # plt.plot(x1,z1)
            # plt.plot(x,z)
            # plt.show()
            if self.spiral:
                if self.centre < self.w0:
                    x = -x
                elif self.centre > self.w0:
                    pass
                elif np.random.rand() > 0.5:
                    x = -x
            elif np.random.rand() > 0.5:
                x = -x
            x += x0
            z += z0
            x = x.astype('int')
            z = z.astype('int')
            xr = np.arange(self.L)
            yr = [np.sum(self.grains[i]) for i in xr]
            for k in xr:
                if k == 0:
                    pass
                else:
                    if x[k] > self.L - 2 or x[k] < 1:
                        self.removed.append(gsize)
                        return 0
                    elif z[k] < yr[x[k]]:
                        # print x[k], v0[layer]
                        # self.removed.append(gsize)
                        self.grains[x[k]].append(gsize)
                        self.z[:,x[k]] += 1 * gsize
                        self.z[:,x[k]-1][1] -= 1 * gsize
                        self.z[:,x[k]+1][0] -= 1 * gsize
                        self.zc[x[k]] =(2*self.grains[x[k]][-2]) - (1*self.grains[x[k]][-1]) + np.random.randint(2,4)
                        return 0
    # @profile
    def drag_force(self,v,a):
        c = self.cd(v,a)
        rho_air = 1.2041
        f = 0.5 * rho_air * np.pi * (a**2) * (v**2) * c
        return f

    def cd(self,v,a):
        vis = 1.51e-5
        Re = float(2. * a * v)/vis
        cd_raw = ((24./Re) * (1. + (0.27 * Re)) **  0.43) + (0.47 * (1. - np.exp(-0.04 * Re ** 0.38)))
        cd = 0.5 + (2.* cd_raw)
        return cd

    def tra(self,gsize,v0,theta):
        theta /= 180.
        theta *= np.pi
        if gsize == 1:
            a = 0.0001
        else:
            a = 0.00075
        n = 400
        dt = 0.0005
        rho = 2650.
        m = (4./3.) * np.pi * (a ** 3) * rho
        x = np.zeros(n,dtype='float')
        z = np.zeros(n,dtype='float')
        v = np.zeros(n,dtype='float')
        vx = np.zeros(n,dtype='float')
        vz = np.zeros(n,dtype='float')
        v[0] = v0
        vx[0] = v0 * np.cos(theta)
        vz[0] = v0 * np.sin(theta)
        for i in range(1,n):
            dvx = self.drag_force(v[i-1],a) * np.absolute(np.cos(theta)) * dt * np.sign(vx[i-1]) / m
            dvz = self.drag_force(v[i-1],a) * np.absolute(np.sin(theta)) * dt * np.sign(vz[i-1]) / m
            # print vx[i-1], dvx, theta
            # print vz[i-1],dvz
            vx[i] = vx[i-1] - dvx
            vz[i] = vz[i-1] - (9.81 * dt) - dvz
            # print (9.81 * dt) + dvz
            v[i] = ((float(vx[i]) ** 2.) + (float(vz[i]) ** 2)) ** 0.5
            theta = np.arctan(vz[i]/vx[i])
            x[i] = x[i-1] + vx[i] * dt
            z[i] = z[i-1] + vz[i] * dt
        x *= 10
        z *= 40
        return x*100,z*100

    def traj(self,size = 1,v0=100.,theta=50.):
        n = 101
        g = 981
        theta *= (np.pi/180.)
        vx = v0 * np.cos(theta)
        vz = v0 * np.sin(theta)
        t = np.linspace(0,0.2,101)
        if size == 1:
            vt = 150.
        else:
            vt = 1000.
        if self.grain_throw == 'drag':
            x = ((v0 * vt * np.cos(theta))/g)*(1 - np.exp(-((g * t)/vt)))
            z = (vt/g)*(v0 * np.sin(theta) + vt)*(1 - np.exp(-((g * t)/vt))) - (vt * t)
        else:
            x = v0 * np.cos(theta) * t
            z = (v0 * np.sin(theta) * t) - ((g/2)*(t**2))
        x *= 10
        z *= 40
        return x,z

# a = Oslo(200,w = 5,p_large = 0.1,grain_throw = 'drag')
# a.run(200)
