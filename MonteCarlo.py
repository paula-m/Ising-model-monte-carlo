""" The checkpoint 1 made by Paulina Mrozek
Date: 21/01/2022

will require boundary conditions
"""

import numpy as np
import random


class Dynamics():
    def __init__(self, N):
        self.N = N

    def energy_from_spin(self, lattice, i, j):
        """instance method calculating the energy of the 4 neighbours from the class."""
        #energy_sum = lattice[i,j]*lattice[i, (j+1)%self.N] + lattice[i, j]*lattice[i, (j-1)%self.N] + lattice[i, j]*lattice[(i+1)%self.N, j] + lattice[i, j]*lattice[(i-1)%self.N, j]
        energy_sum = lattice[i,j]*(lattice[i, (j+1)%self.N] + lattice[i, (j-1)%self.N] + lattice[(i+1)%self.N, j] + lattice[(i-1)%self.N, j])   #the energy created
        return -energy_sum

    def glauber_update(self, m, spin_lattice, T):
        """Method for update of glauber dynamics - 1 metropolis sweep"""
        for n in range(m):  #m is number of iterations for 1 full sweep
            rsi, rsj = random.randrange(self.N), random.randrange(self.N)   #selects random sites by selecting values in lattice size
            delta_E = -2*self.energy_from_spin(spin_lattice, rsi, rsj)  #calculates energy change for the flip
            if delta_E <= 0:    #acceptance condition
                spin_lattice[rsi, rsj] = -spin_lattice[rsi, rsj]
            else:
                if random.uniform(0, 1) <= np.exp(-delta_E/T):  #if not initial accepted then accept with probability min(1, -deltaE/t)
                    spin_lattice[rsi, rsj] = -spin_lattice[rsi, rsj]
        return spin_lattice     #returns final lattice after the sweep.

    def kawasaki_energy_next(self, spin_lattice, x1, y1, x2, y2, T):
        """Method for when the sites are next to each other"""
        totE = -2*self.energy_from_spin(spin_lattice, x2, y2) + -2*self.energy_from_spin(spin_lattice, x1, y1) + 4  #accounts for the overcouting
        if totE <=0:
            spin_lattice[x1, y1] = -spin_lattice[x1, y1]
            spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        else:
            if random.uniform(0, 1) <= np.exp(-totE/T):     #changes the sign of both of them ===> switches the two spins. 
                spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        return spin_lattice


    def kawasaki_normal_energy_change(self, spin_lattice, x1, y1, x2, y2, T):
        """Method for normal kawasaki energy change"""
        totE = -2*self.energy_from_spin(spin_lattice, x2, y2) + -2*self.energy_from_spin(spin_lattice, x1, y1)  #counting for the two changes for kawasaki
        if totE <= 0:
            spin_lattice[x1, y1] = -spin_lattice[x1, y1]
            spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        else:
            if random.uniform(0, 1) <= np.exp(-totE/T):
                spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        return spin_lattice



    def kawasaki_dynamics(self, n, m, spin_lattice, T):
        """method for updating the kawasaki dynamics"""
        i = 0   #the counting of the metropolis
        while i <= m:   #making sure that we are only counting if the are not the same spin hence why while and not if or for loops
            x1, x2, y1, y2 = random.randrange(n), random.randrange(n), random.randrange(n), random.randrange(n) #creates the random x, y vals for the randomly selected spins.
            if spin_lattice[x1, y1] != spin_lattice[x2, y2]:    #checks if they are the same spin
                i += 1  #spin is counted only if they are different spins.
                if ((abs(x2 - x1) == 1) & (y1 == y2)) | ((abs(y1 - y2) ==1) & (x1 == x2)):  #checks if they are next to each other
                    spin_lattice = self.kawasaki_energy_next(spin_lattice, x1, y1, x2, y2, T)   #if yes then account for it
                else:
                    spin_lattice = self.kawasaki_normal_energy_change(spin_lattice, x1, y1, x2, y2, T)  #if not then normal energy change
        return spin_lattice
