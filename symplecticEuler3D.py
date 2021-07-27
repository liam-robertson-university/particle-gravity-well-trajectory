"""
CMod Ex3: symplectic Euler time integration of
particle in gravitational potential.

Produce a plot of the particles movement in the x-y plane and an output file containing the
x-y-z coordinates of the particle with time.

The format of the output file is (by line): xpos ypos zpos time.
"""

import sys
import numpy as np
import math
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

# Read name of output file from command line
if len(sys.argv)!=2:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + " <output file>"
    quit()
else:
    outfileName = sys.argv[1]

# Open output file for writing
outfile = open(outfileName, "w")

# Define a function to calculate the gravitational force between 2 particles
# input: mass_c is the mass of a particle at the origin (float), p is a Particle3D instance
# output: gravitational force acting on p (an array)
def calculateGForce(mass_c, p):
    return -p.mass*mass_c*p.position/(p.magPos()**3)

# Define a function to calculate the gravitational potential between 2 particles
# input: mass_c is the mass of a particle at the origin (float), p is a Particle3D instance
# output: gravitational potential of p (float)
def calculateGPot(mass_c, p):
    return -p.mass*mass_c/p.magPos()

# Set up particle
mass_p = 0.0001
pos  = np.array([1,0,0],float)
vel  = np.array([0,1,0],float)
label = "particle in gravity"
p = Particle3D(pos, vel, mass_p, label)

# Set up simulation parameters
numstep = 500
time = 0.
dt = 0.1
mass_c = 1.

# Set up initial total energy
initE = p.kineticEnergy() + calculateGPot(mass_c, p)

# Set up data lists
xValue = [p.position[0]]
yValue = [p.position[1]]
zValue = [p.position[2]]
tValue = [time]
deltaElist = [] # This list is used to find the highest energy fluctuation in the simulation

# Write initial particle data to the file
outfile.write("{0:f} {1:f} {2:f} {3:f}\n".format(p.position[0], p.position[1], p.position[2], time))

# Start the time integration loop

for i in range(numstep):
    # Update particle position
    p.leapPos1st(dt)
    # Update force
    force = calculateGForce(mass_c, p)
    # Update particle velocity
    p.leapVelocity(dt, force)

    # Increase time
    time = time + dt

    # Update total energy
    totalE = p.kineticEnergy() + calculateGPot(mass_c, p)

    # Append (magnitude of) current energy fluctuation to list
    deltaElist.append(((initE - totalE)**2)**(1.0/2.0))
    
    # Append particle data to the relevant lists
    xValue.append(p.position[0])
    yValue.append(p.position[1])
    zValue.append(p.position[2])
    tValue.append(time)
    outfile.write("{0:f} {1:f} {2:f} {3:f}\n".format(p.position[0], p.position[1], p.position[2], time))

# Find the highest energy fluctuation using a list command and display it
maxdeltaE = max(deltaElist)
print(maxdeltaE)

# Close output file
outfile.close()

# Plot graph of the particles position projected in the x-y plane
pyplot.plot(xValue,yValue)
pyplot.title("x-y position of a particle under gravity")
pyplot.xlabel("x-position (m)")
pyplot.ylabel("y-position (m)")
pyplot.show()
