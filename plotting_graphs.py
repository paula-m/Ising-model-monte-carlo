from turtle import color
import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

class plotting_graphs():

    def plot_glauber(self):
        column_names = ["Temperature", "Magnetisation", "Energy", "Susceptibility", "Specific", "Errorheat", "Errorsus"]
        df = pd.read_csv("data_glauber_t1.csv", names = column_names)
        Energy = df.Energy.to_list()
        Temp = df.Temperature.to_list()
        Magnetisation = df.Magnetisation.to_list()

        Suscept = df.Susceptibility.to_list()
        Heat_cap = df.Specific.to_list()
        Heaterr = df.Errorheat.to_list()
        Suserr = df.Errorsus.to_list()
        

        plt.plot(Temp[0:21], Energy[0:21], label="1st data")
        plt.title("The Energy vs Temp - Glauber")
        plt.xlabel("Temperature")
        plt.ylabel("Energy")
        plt.legend()
        plt.show()

        plt.plot(Temp[0:21], Magnetisation[0:21], label="Normal magnetisation")
        plt.title("Magnetisation vs Temp - Glauber")
        plt.xlabel("Temperature")
        plt.ylabel("Magnetisation")
        plt.legend()
        plt.show()
        
        plt.plot(Temp[0:21], Suscept[0:21], color="orange", ms=5)
        plt.errorbar(Temp[0:21], Suscept[0:21], yerr=Suserr[0:21], ecolor="purple", fmt="none")
        plt.title("Susceptibility vs Temp - Glauber")
        plt.xlabel("Temperature")
        plt.ylabel("Susceptibility")
        plt.show()
        
        plt.plot(Temp[0:21], Heat_cap[0:21])
        plt.errorbar(Temp[0:21], Heat_cap[0:21], yerr=Heaterr[0:21], ecolor="purple", fmt="none")
        plt.title("Heat Capacity vs Temperature - Glauber")
        plt.xlabel("Temperature")
        plt.ylabel("Heat Capacity")
        plt.show()
        

    def plot_kawasaki(self):
        column_names = ["Temperature", "Energy", "Specific", "Error"]
        df = pd.read_csv("data_kawasaki.csv", names=column_names)

        
        energy = df.Energy.to_list()
        specific = df.Specific.to_list()
        Temperature = df.Temperature.to_list()
        Error = df.Error.to_list()


        plt.plot(Temperature, energy)
        plt.title("the temperature vs energy - Kawasaki")
        plt.xlabel("Temperature")
        plt.ylabel("Energy")
        plt.show()
        
        plt.plot(Temperature, specific)
        plt.errorbar(Temperature, specific, yerr=Error)
        plt.title("The Heat Capacity vs Temp - Kawasaki")
        plt.xlabel("Temperature")
        plt.ylabel("Heat Capacity")
        plt.show()
    
def main():
    a = plotting_graphs()
    a.plot_kawasaki()
    a.plot_glauber()
main()