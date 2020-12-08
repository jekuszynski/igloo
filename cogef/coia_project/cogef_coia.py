from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.calculators.morse import MorsePotential #change this to whatever calculator you want to use
from cogef import COGEF, Dissociation
import pylab as plt #cute shortcut for matplotlib~ (and numpy?)

def parse_gjf(file_name,atom_count):
    f=open(file_name) #open file
    lines = f.readlines() #read file contents as 'lines'

    l = [None]*len(lines) #make empty list with # of total lines
    l2 = [None]*atom_count #make empty list with number of atoms
    j=0 #start with index = 0
    for line in lines: #repeat below for every # of lines
        l[j] = line.split() #split line at every index
        j+=1 #add to index by 1

    elements=[] #make empty list for atoms
    xyz=[] #make empty list for coordinates
    for i in range(len(l2)): #for every line which has an atom
        l2[i] = l[i+6] #skip x amount of lines until atoms are reached in .gjf file
        elements.append(l2[i][0]) #add element to elements list
        x = float(l2[i][1]) #add x to xyz
        y = float(l2[i][2]) #add y to xyz
        z = float(l2[i][3]) #add z to xyz
        xyz.append((x,y,z)) #add all xyz to xyz list
    f.close() #close file

    return elements,xyz #return function with elements and xyz lists

def initialize(image, imagenum, new_opt, get_filename): #function required for cogef.pull 
    """Initialize the image and return the trajectory name."""
    if get_filename:
        return 'cogef_coia' + str(imagenum) + '.traj'
    image.set_calculator(MorsePotential())

#---Executable Code Below---#

'''Gives element list and positions list from .gjf file
Second argument (e.g. '42'), should be entered for number of atoms in .gjf file'''
elements,positions = parse_gjf("coGEF1.gjf",42)

'''Prepare 'image' for cogef using calculator and FIRE'''
image = Atoms(elements, positions)
image.set_calculator(MorsePotential()) #set calculation theory. Morse used as placeholder.
fmax = 0.05 
FIRE(image).run(fmax=fmax)

'''Relax atoms using cogef'''
cogef = COGEF([image], 34, 38, optimizer=FIRE, fmax=fmax)

'''Pull apart atoms using cogef'''
stepsize = 0.05
steps = 30
cogef.pull(stepsize, steps, initialize, 'cogef_coia.traj')

'''Set important data to variables and plot'''
energies, distances = cogef.get_energy_curve() #calculate energy and distance from cogef run
plt.figure(0) #prepare figure for plotting
plt.xlabel('Displacement from Equilibrium, D ($\\AA$)') #set x axis label
plt.ylabel('U [eV]') #set y axis label
plt.plot(distances,energies) #plot distance v. energy
plt.savefig('coia_data.png',dpi=300) #can also change to path (e.g. 'Sub Directory/graph.png') to set location for saved plot)