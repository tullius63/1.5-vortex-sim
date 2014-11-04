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

for i in range(0,N):
    posx[i] = L*random()
    posy[i] = L*random()
 
for i in range(0,5000):
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
    energies.append(vsum)          
					# Update regardless of acceptance!

#r2 = n.zeros(N-1)
#arrcount = n.zeros((N,3588)) 
#avecount = n.zeros(3588)
#for i in range(0,N):
#    print i
#    p0x = posx[i]
#    p0y = posy[i]
#    cunt = -1
#    count = 0
#    for j in range(0,N):
#        if j is not i:  
#            cunt+=1
#            r2[cunt] = m.sqrt((posx[cunt]-p0x)**2+(posy[cunt]-p0y)**2)
#        else:
#            continue
#    for k in range(0,3588):
#        for l in range(0,N-1):
#            if r2[l] < (.125+ k*.01):
#                if r2[l] >= (.125+ (k-1)*.01):
#                    count+=1
#        arrcount[l,k]=count
#for i in range(0,3588):
#    rdf[i] = n.mean(arrcount[:,i])
#xaxis= list(n.arange(.125, 36.05, .01))




 #Calculating G(r) (radial distribution function) of a water box, 
 #taking into account periodic boundaries
histo = n.zeros((N,3000))
r2 = n.zeros(N-1)
r =list(range(0,3000))
pos = n.vstack([posx, posy]).T
pos_periodic = n.vstack(pos+displ for displ in [[-L,-L], [-L,0], 
                                                 [-L,+L], [+L,-L],
                                                 [+L,0], [+L, +L],
                                                 [0,-L], [0, +L],
                                                 [0,0]])
                                                  
rad = scipy.spatial.distance.pdist(pos_periodic)



rad1 = scipy.spatial.distance.squareform(rad)



rad1 = n.ma.array(rad1)
rad1 = n.ma.masked_equal(rad1,0)

print rad1
nearest = n.amin(rad1, axis=1)

#xaxis= list(n.arange(rdf.min(), rdf.max(), .011958))
#fig2 = plt.figure(1)
#fig2 = plt.scatter(xaxis, rdf)
#fig2 = plt.xlabel('r', fontsize=18)
#fig2 = plt.ylabel('g(r)', fontsize=18)

data, bins = n.histogram(rad, bins = 3000)

dr = bins[1] - bins[0]
g = data/dr/2/n.pi/bins[:-1]/N
            
            

#ExportName_csv = 'C:\Users\Sam\Desktop\ExportData/' + N + T + '.csv'
#data = zip(*np.histogram(rad,3000))
#np.savetxt(ExportName_csv, izip(freq, bins), delimiter=",")        
    

#Normalize RDF

fig1=plt.figure(1)
plt.clf()            
plt.scatter(posx,posy)
plt.xlabel('x')
plt.ylabel('y')


fig2 = plt.figure(2)
fig2.clf()   
plt.plot(bins[:-1], g)
plt.xlim(0,L)
plt.xlabel('r', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
   

fig3 = plt.figure(3)
fig3.clf()   
fig3 = plt.plot(energies)

fig4 = plt.figure(4)
fig4.clf()   
fig4 = plt.hist(nearest, bins = 3000, normed = True)

plt.show() 