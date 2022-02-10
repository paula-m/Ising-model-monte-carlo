""" The checkpoint 1 made by Paulina Mrozek
Date: 21/01/2022

will require boundary conditions
"""
from turtle import color
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
matplotlib.use('TKAgg')

class Glauber_dynamics():
    def __init__(self, N):
        self.N = N

    def creating_array(self):
        spin_lattice = np.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                spin_lattice[i,j] = random.choice([-1,1])
        print(spin_lattice.shape)
        return spin_lattice


    def change_in_energy(self, lattice, i, j):
        energy_sum = lattice[i,j]*lattice[i, (j+1)%self.N] + lattice[i, j]*lattice[i, (j-1)%self.N] + lattice[i, j]*lattice[(i+1)%self.N, j] + lattice[i, j]*lattice[(i-1)%self.N, j]
        return -energy_sum

    def iterate_over_array(self, m, spin_lattice):
        print(np.count_nonzero(spin_lattice == 1))
        for n in range(m):
            rsi, rsj = random.randrange(self.N), random.randrange(self.N)
            #energy_init = self.change_in_energy(spin_lattice, rsi, rsj)
            delta_E = -2*self.change_in_energy(spin_lattice, rsi, rsj)
            if delta_E <= 0:
                spin_lattice[rsi, rsj] = -spin_lattice[rsi, rsj]
            else:
                if random.uniform(0, 1) <= np.exp(-delta_E):
                    spin_lattice[rsi, rsj] = -spin_lattice[rsi, rsj]
        return spin_lattice

    def kawasaki_energy(self, spin_lattice, x1, y1, x2, y2):
        energy1 = self.change_in_energy(spin_lattice, x1, y1)
        delta_e1 = -2*energy1
        energy2 = self.change_in_energy(spin_lattice, x2, y2)
        delta_e2 = -2*energy2
        totE = -2*energy2 + delta_e1 + 2
        if totE <=0:
            spin_lattice[x1, y1] = -spin_lattice[x1, y1]
            spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        else:
            u = random.uniform(0, 1)
            if u <= np.exp(-totE):
                spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                spin_lattice[x2, y2] = -spin_lattice[x2, y2]
        return spin_lattice

    def kawasaki_dynamics(self, n, m, spin_lattice):
        print( np.count_nonzero(spin_lattice == 1))

        im = plt.imshow(spin_lattice, animated = True, cmap='Pastel1')


        #for i in range(m):
        i = 0
        while i <= m:
            x1, x2, y1, y2 = random.randrange(n), random.randrange(n), random.randrange(n), random.randrange(n)

            #spin1 = spin_lattice[x1, y1]
            #spin2 = spin_lattice[x2, y2]
            print("x1a", x1, "x2a", x2, "y1", y1, "y2", y2)
    

            if spin_lattice[x1, y1] != spin_lattice[x2, y2]:
                i += 1
                if ((x1 == (x2 +1)%n) & (y1 == y2)) | ((x1 == (x2-1)%n) & (y1 == y2)) | ((y1 == (y2-1)%n) & (x1 == x2)) | ((y1 == (y2 +1)%n) & (x1 == x2)):
                    print("x1", x1, "y1", y1, "x2", x2, "y2", y2)
                    spin_lattice = self.kawasaki_energy(spin_lattice, x1, y1, x2, y2)
                    #energy1 = self.change_in_energy(spin_lattice, x1, y1)
                    #delta_e1 = -2*energy1
                    #energy2 = self.change_in_energy(spin_lattice, x2, y2)
                    #delta_e2 = -2*energy2
                    #totE = -2*energy2 + delta_e1 + 2
                    #if totE <=0:
                    #    spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                    #    spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                    #else:
                    #    u = random.uniform(0, 1)
                    #    if u <= np.exp(-totE/5):
                    #        spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                    #        spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                    #    else:
                    #        continue
                #elif ((y1 == (y2 +1)%n) & (x1 == x2)) | ((y1 == (y2-1)%n) & (x1 == x2)):
                #    energy1 = self.change_in_energy(spin_lattice, x1, y1)
                #    delta_e1 = -2*energy1
                #    energy2 = self.change_in_energy(spin_lattice, x2, y2)
                #    delta_e2 = -2*energy2
                #    totE = delta_e2 + delta_e1 +2
                #    if totE <=0:
                #        spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                #        spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                #    else:
                #        u = random.uniform(0, 1)
                #        if u <= np.exp(-totE/5):
                #            spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                #            spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                #        else:
                #            continue
                else:
                    spin_lattice = self.kawasaki_energy(spin_lattice, x1, y1, x2, y2)
                    #energy1 = self.change_in_energy(spin_lattice, x1, y1)
                    #delta_e1 = -2*energy1
                    #energy2 = self.change_in_energy(spin_lattice, x2, y2)
                    #delta_e2 = -2*energy2
                    #totE = delta_e2 + delta_e1
                    #if totE <=0:
                    #    spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                    #    spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                    #else:
                    #    u = random.uniform(0, 1)
                    #    if u <= np.exp(-totE/5):
                    #        spin_lattice[x1, y1] = -spin_lattice[x1, y1]
                    #        spin_lattice[x2, y2] = -spin_lattice[x2, y2]
                    #    else:
                    #        continue
            #print(i)
            if i%10 == 0:
                f = open('spins.data', 'w')
                for k in range(len(spin_lattice[0])):
                    for l in range(len(spin_lattice)):
                        f.write('%d %d %lf \n'%(k,l, spin_lattice[k,l]))
                f.close()
                im.set_data(spin_lattice)
                plt.draw()
                plt.pause(0.001)
        print("mean", np.mean(spin_lattice))
        return spin_lattice


        
        
        #return 

    def saving_data(self, nosteps, m, spin_lattice):
        print(spin_lattice)
        #fig = plt.figure()
        im = plt.imshow(spin_lattice, animated = True, cmap='Pastel1')
        for i in range(nosteps):
            #spin_lattice = self.iterate_over_array(m, spin_lattice)
            spin_lattice = self.kawasaki_dynamics(len(spin_lattice), m, spin_lattice)
            #f = open('spins.data', 'w')
            #for k in range(len(spin_lattice[0])):
            #    for l in range(len(spin_lattice)):
            #        f.write('%d %d %lf \n'%(k,l, spin_lattice[k,l]))
            #f.close()
            #im.set_data(spin_lattice)
            #plt.draw()
            #plt.pause(0.001)
        #plt.close()
        return spin_lattice
        

def main():
    a = Glauber_dynamics(50)
    spin_lattice = a.creating_array()
    #print_the_thing = a.iterate_over_array(1000, spin_lattice)
    #saving = a.saving_data(1000, 2500, spin_lattice)
    #kawa  = a.kawasaki_dynamics(50, 100, spin_lattice)
    anim = a.saving_data(2500, 2500, spin_lattice)


main()
