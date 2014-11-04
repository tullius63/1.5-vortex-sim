import math as m
from random import *
import numpy as n
import matplotlib.pyplot as plt
import scipy.spatial.distance
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


N= 300
L = 36
T = 1
beta = 1/(1.3806488*10**-23*T)
posx = n.zeros(N)
posy = n.zeros(N)
energies = []

def pot(r):                                #Potential function
    pot = 11.33*(m.exp(-0.8645*r**2) - 0.9846*m.exp(-0.8534*r**2))/(11.33*(1-.9846))
    return pot

def simulate(posx0, posy0, beta, N=300, steps=20000):
    posx, posy = posx0, posy0
    for i in range(0,20000):
        vsum = 0
        index= randint(0,N-1)
        rx = posx[index]
        ry = posy[index]
        v=0
        vnew=0
    
        theta= 2*m.pi*random()
        step = .125 + 3*random()/8
        rx_new =rx + step*m.cos(theta)
        ry_new =ry + step*m.sin(theta)
        if rx_new > L:
            rx_new= rx_new -L
        if rx_new < 0:
            rx_new = L + rx_new
        if ry_new > L:
            ry_new = ry_new -L
        if ry_new < 0:
            ry_new = L - ry_new
        for i in range(0,N):
            if i is index:
                continue
            else:
                r = m.sqrt((rx-posx[i])**2 + (ry-posy[i])**2)    
                v=v+pot(r)
                r_new= m.sqrt((rx_new-posx[i])**2 + (ry_new-posy[i])**2)    
                vnew=vnew+pot(r_new)
                
        if vnew < v:				# Downhill move always accepted
            posx[index] = rx_new
            posy[index] = ry_new
            v = vnew
        else:					# Uphill move accepted with luck
            A_mov = m.exp( -beta *(vnew - v) )
            if A_mov > random():
                posx[index] = rx_new
                posy[index] = ry_new
                v = vnew
        vsum = vsum + v
        energies.append(vsum)          # Update regardless of acceptance!
        
    pos = n.vstack([posx, posy]).T
    return pos, energies

def run_simulation(beta):
    pos0 = L * n.random.random((N, 2))
    pos, _ = simulate(pos0[:,0], pos0[:,1], beta)
    return pos
    
def generate(n=10):
    configs = {}
    import multiprocessing
    pool = multiprocessing.Pool(4)
    configs = pool.map(run_simulation, range(4))
    configs = dict(enumerate(configs))
    n.io.savez('configs.npz', **configs)

def analyze(configs):
    r_range = (0, L)
    rdfs = []
    for pos in configs:
        #Calculating G(r) (radial distribution function) of a water box, 
        #taking into account periodic boundaries
        pos_periodic = n.vstack(pos+displ for displ in [[-L,-L], [-L,0], 
                                                        [-L,+L], [+L,-L],
                                                        [+L,0], [+L, +L],
                                                        [0,-L], [0, +L],
                                                        [0,0]])                                                    
        rad = scipy.spatial.distance.pdist(pos_periodic)
        
        #xaxis= list(n.arange(rdf.min(), rdf.max(), .011958))
        #fig2 = plt.figure(1)
        #fig2 = plt.scatter(xaxis, rdf)
        #fig2 = plt.xlabel('r', fontsize=18)
        #fig2 = plt.ylabel('g(r)', fontsize=18)
        
        data, bins = n.histogram(rad, bins = 3000, range=r_range)
        
        dr = bins[1] - bins[0]
        g = data/dr/2/n.pi/bins[:-1]/N
        rdfs.append((bins, g))
        
        #ExportName_csv = 'C:\Users\Sam\Desktop\ExportData/' + N + T + '.csv'
        #data = zip(*np.histogram(rad,3000))
        #np.savetxt(ExportName_csv, izip(freq, bins), delimiter=",")        
        
        #plt.scatter(posx,posy)
        #plt.xlabel('x')
        #plt.ylabel('y')
    
    rdf = plt.figure(1).gca()
    rdf.cla()
    for bins, g in rdfs:
        rdf.plot(bins[:-1], g, '+')
        
    rdf.errbar(bins[:-1], y=np.mean(rdfs, axis=1), yerr=np.std(rdfs, axis=1))
    rdf.set_xlabel('r', fontsize=18)
    rdf.set_ylabel('g(r)', fontsize=18)
    rdf.show() 

if __name__ == '__main__':    
    generate()
    #configs = n.io.load('configs.npz')
    #analyze(configs)