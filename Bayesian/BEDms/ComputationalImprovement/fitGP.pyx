import numpy as np

import matplotlib.pylab as plt

import gaussianprocess as gp



x= np.linspace(0, 1, 20)



y= []

with open('PartInferenceResults_DynStim_11_Simulation_DynStim_11_GFP.csv') as f:

    for line in f:

        if 'V' not in line and len(line) > 1:

            y.append(line.rstrip().split(',')[1:])

y= np.array(y).astype(float)

t= np.arange(y.shape[0])



# convert to 1D arrays

x= np.tile(t, y.shape[1])

y= y.flatten(order= 'F')



# reduce number of data points

n= 500

kp= np.random.randint(len(y), size= n)

xs, ys= x[kp], y[kp]

ys= ys[np.argsort(xs)]

xs= np.sort(xs)



# fit GP

g= gp.maternGP({0: (0, 8), 1: (-4, 4), 2: (0, 6)}, xs, ys)

g.findhyperparameters(noruns= 5)

g.results()

g.predict(xs)

plt.figure()

g.sketch('.')

plt.plot(xs, g.sample(5), '.-')

plt.show()



# find entropy

gcv= 2*np.pi*np.e*g.covp

# see pseudo-determinant in Wikipedia

from scipy.linalg import svd

s= svd(gcv, compute_uv= False)

ent= np.sum(np.log(s))/2

print(ent)

