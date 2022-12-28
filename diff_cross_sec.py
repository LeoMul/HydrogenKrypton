import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy import special as s
import cmath

data = np.loadtxt("DataFiles/phasesEmin0.100Emax3.200N009996h0.000010.dat")
energies = data[:,0]
l_array = (0,1,2,3,4,5,6)
theta = np.linspace(0,np.pi,100)
costheta = np.cos(theta)
mu = 1.6537e-27
hbar = 1.05457182e-34
#columns are p_l(costheta)
legendre_matrix = np.zeros([len(costheta),len(l_array)])

for l in range(0,len(l_array)):
    x = s.legendre(l)
    legendre_matrix[:,l] = x(costheta)

#print(energies)

def make_cross_sec_matrix():
    dcsa_matrix = np.zeros([len(theta),len(energies)])

    for i in range(0,len(energies)):
        kvectorsquared = energies[i] * 1.60218E-22 *1e-20* 2*mu/(hbar**2) 
        for p in range(0,len(theta)):
            summation = 0.0
            for l in range(0,len(l_array)):
                phase = data[i,l+1] 
                summation = summation + (2*l_array[l]+1)* np.sin(phase)*legendre_matrix[p,l]*cmath.exp(1.0j*phase)
            summation = np.abs(summation)**2 
            dcsa_matrix[p,i] = summation/kvectorsquared
    return dcsa_matrix

dcsa_matrix = make_cross_sec_matrix()
print(dcsa_matrix)

def solid_angle_integrate(diffcsa,theta):
    integral = 0
    array = np.sin(theta)
    array = array * diffcsa 
    dtheta = theta[1]-theta[0]
    integral = sum(array)*dtheta * np.pi*2

    return integral

def create_csa():
     cross_sec_vec = np.zeros(len(energies))
     for j in range(0,len(cross_sec_vec)):
         dcsa = dcsa_matrix[:,j]
         cross_sec_vec[j] = solid_angle_integrate(dcsa,theta)

     return cross_sec_vec

#csa_vec = create_csa()

#plt.figure()
#plt.plot(energies,csa_vec)
#plt.show()

def create_plot():
        plt.rcParams['text.usetex'] = True
        T_array, X_array = np.meshgrid(energies, theta  )
        x = [0,np.pi/4,np.pi/2,0.75*np.pi,np.pi]
        labels = ['0','${\pi}/{4}$','$\pi/2$','${3\pi}/{4}$','$\pi$']
        y = [0,0.875, 1.75,2.625, 3.5]
        plt.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(X_array, T_array, dcsa_matrix, cmap=cm.coolwarm,linewidth=0, antialiased=False)
        ax.zaxis.set_major_locator(LinearLocator(10))
        plt.xlabel(r"$\theta$/rad", labelpad=10)
        plt.ylabel("$E$/meV", labelpad=10)
        ax.zaxis.set_rotate_label(False) 
        ax.view_init(elev=17., azim = -56)
        ax.set_zlabel(r"$d\sigma/d\Omega$",rotation = 0.0, labelpad=10)
        plt.xticks(x,labels)
        ax.set_yticks(y)
        ax.set_zticks([0, 250,500, 750,1000])
        fig.colorbar(surf, shrink=0.5, aspect=5,orientation="vertical", pad = 0.18)
        plt.savefig("diffcrosssec.pdf")
        plt.show()




create_plot()