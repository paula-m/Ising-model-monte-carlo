from ensurepip import bootstrap
from MonteCarlo import Dynamics

import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
matplotlib.use('TKAgg')
import csv


class MonteCarlo():
    """The Class which goes through the metropolis algorithm, visualisation and calculation of the errors"""
    def __init__(self):
        self.N = int(input("size "))
        self.dynamics = Dynamics(self.N) #initialises the class for the dynamics
    
    def creating_array(self):
        """The method which creates random array for the random start."""
        spin_lattice = np.zeros((self.N,self.N))    #creates numpy array of zeros
        for i in range(self.N):
            for j in range(self.N):
                spin_lattice[i,j] = random.choice([-1,1])   #replacing each number with random number 1 or -1
        return spin_lattice
    
    
    def animation(self, spin_lattice):
        f = open('spins.data', 'w') #opens the data for animation
        for k in range(len(spin_lattice[0])):
            for l in range(len(spin_lattice)):
                f.write('%d %d %lf \n'%(k,l, spin_lattice[k,l]))    #writesthe spin and coordinates
        f.close()   #closes file
        plt.cla()   #clears the figure plot to plot the next animation
        im=plt.imshow(spin_lattice, animated=True, cmap='Pastel1')  #animates the spin lattice
        plt.draw()
        plt.pause(0.001)

        
    def monte_carlo_update_kawasaki(self, nosteps, m, spin_lattice, T):
        """The monte carlo kawasaki dynamics - we are not looking at magnetisation and susceptibility because it stays the same."""
        
        enery_arr = []
        
        for i in range(nosteps):    #the sweep (MC)
            spin_lattice = self.dynamics.kawasaki_dynamics(len(spin_lattice), m, spin_lattice, T)   #calling the kawasaki dynamics
            
            self.animation(spin_lattice)    #calls the animation function comment out if animation not needed
        
            if i >= 200:    #looking at the metropolis sweeps after they already have been thermalised
                if (i%10) == 0: #only want to get measurements when they are independent each 10 sweeps
                    energy = 0  #calculates energy
                    for p in range(len(spin_lattice[0])):
                        for n in range(len(spin_lattice)):
                            energy += self.dynamics.energy_from_spin(spin_lattice, p, n)    #this calculates energy from each
                    enery_arr.append((energy/2))    #to account for overcounting and boundary conditions.

        scaled_heat_cap = (1/(m*(T**2)))*((np.mean(np.array(enery_arr)**2))-((np.mean(enery_arr))**2))

        bootstrap_heat = self.bootstrap_specific_heat(enery_arr, T, m) #calls the function to get the error on the heat capacity for this temp
        
        k = open('data_kawasaki.csv', 'a+')  #saves the data from kawasaki for a csv file for calculations later on.
        writer = csv.writer(k)
        writer.writerow([T, np.mean(enery_arr), scaled_heat_cap, bootstrap_heat])
        k.close()

        return spin_lattice
        
        
    def monte_carlo_update_glauber(self, nosteps, m, spin_lattice, T):
        """method for the glauber update same type as kawasaki but has the magnetisation"""
        magnetisation_array = []
        enery_arr = []

        for i in range(nosteps):
            spin_lattice = self.dynamics.glauber_update(m, spin_lattice, T)
            self.animation(spin_lattice)
            
            if i >= 200:
                if (i%10) == 0:
                    a = np.sum(spin_lattice)
                    magnetisation_array.append(a)
                    energy = 0
                    for p in range(len(spin_lattice[0])):
                        for n in range(len(spin_lattice)):
                            energy += self.dynamics.energy_from_spin(spin_lattice, p, n)
                    enery_arr.append((energy/2))
        plt.close()
        susceptibility = (1/(m*T))*((np.mean(np.array(magnetisation_array)**2)) - (np.mean(magnetisation_array)**2))
        scaled_heat_cap = (1/(m*(T**2)))*((np.mean(np.array(enery_arr)**2))- (np.mean(enery_arr)**2))

        bootstrap_err_heat = self.bootstrap_specific_heat(enery_arr, T, m)
        bootstrap_susc = self.bootstrap_susceptibility(magnetisation_array, T, m)
        
        k = open('data_glauber_t1.csv', 'a')
        writer = csv.writer(k)
        writer.writerow([T, np.mean(magnetisation_array), np.mean(enery_arr), susceptibility, scaled_heat_cap, bootstrap_err_heat, bootstrap_susc])
        k.close()
        return spin_lattice    


    def choice_of_dynamics(self, nosteps, spin_lattice):
        """The method class for the monte carlo update asking for input"""
        dynamics_choice = input("Please give dynamics, kawasaki or glauber: ")
        T = int(input("Temperature "))
        latticesize = (self.N)**2
        if dynamics_choice == "kawasaki":
            self.monte_carlo_update_kawasaki(nosteps, latticesize, spin_lattice, T)
        elif dynamics_choice == "glauber":
            self.monte_carlo_update_glauber(nosteps, latticesize, spin_lattice, T) 
        else:
            print("Please type correct dynamics")   #if typo then try again
            self.choice_of_dynamics(nosteps, latticesize, spin_lattice, T)
      
    
    def bootstrap_specific_heat(self, lattice, T, m):
        """Errors for calculating specific heat """
        all_c = []
        for j in range(1000):
            new_sample = []
            for i in range(len(lattice)):
                ind = random.randrange(len(lattice))    #selects index of random site
                new_sample.append(lattice[ind]) #appends the selected energy into array
            spec_heat = ((np.mean(np.array(new_sample)**2)) - (np.mean(np.array(new_sample)))**2)*(1/((m)*(T**2)))   #calculates SHC for the new array of energies 
            all_c.append(spec_heat)
        error = np.sqrt((np.mean(np.array(all_c)**2)) - (np.mean(np.array(all_c))**2))  #calculates the error.
        return error
    
    def bootstrap_susceptibility(self, lattice, T, m):
        """Errors for calculating susceptibility"""
        all_c = []
        for j in range(1000):
            new_sample = []
            for i in range(len(lattice)):
                ind = random.randrange(len(lattice))
                new_sample.append(lattice[ind])
            spec_heat = ((np.mean(np.array(new_sample)**2)) - (np.mean(np.array(new_sample)))**2)*(1/((m)*T))
            all_c.append(spec_heat)
        error = np.sqrt((np.mean(np.array(all_c)**2)) - (np.mean(np.array(all_c))**2))
        return error
            
      

def main():
    try_dyn = MonteCarlo()
    

    spin_lattice = try_dyn.creating_array()
    try_dyn.choice_of_dynamics(2500, spin_lattice)
    
    #The code for long run for glauber 
    #spin_lattice = (-1)*np.ones((50,50))
    #temp_arr = np.arange(1, 3.1, 0.1)
    #for j in range(len(temp_arr)):
    #    spin_lattice = try_dyn.monte_carlo_update_glauber(10200, 2500, spin_lattice, temp_arr[j])



    #The code for long run for kawasaki
    #spin_lattice = np.ones((50, 50))
    #for j in range(50):
    #    for k in range(25):
    #        spin_lattice[j][k] = -1

    #for i in range(len(temp_arr)):
        #spin_lattice = try_dyn.monte_carlo_update_kawasaki(10200, 2500, spin_lattice, temp_arr[i])
        

main()
