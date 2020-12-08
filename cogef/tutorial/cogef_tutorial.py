from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.calculators.morse import MorsePotential
from ase.units import J, m
from cogef import COGEF, Dissociation
import pylab as plt #cute shortcut for matplotlib~ (and numpy?)

def initialize(image, imagenum, new_opt, get_filename):
    """Initialize the image and return the trajectory name."""
    if get_filename:
        return 'cogef' + str(imagenum) + '.traj'
    image.set_calculator(MorsePotential())

fmax = 0.05

image = Atoms('H10', positions=[(i, 0, 0) for i in range(10)])
image.set_calculator(MorsePotential()) #set calculation theory. Morse used as placeholder.
FIRE(image).run(fmax=fmax)

cogef = COGEF([image], 0, 8, optimizer=FIRE, fmax=fmax) #relax atoms through cogef using optimizer
stepsize = 0.03
steps = 35

cogef.pull(stepsize, steps, initialize, 'cogef.traj') #simulate pulling apart of atoms

energies, distances = cogef.get_energy_curve() #calculate energy and distance from cogef run
plt.figure(0) #create figure (empty)
plt.plot(distances, energies) #make plot
plt.xlabel('d [$\\AA$]') #set x axis label
plt.ylabel('U [eV]') #set y axis label
plt.savefig('cogef_tutorial.png') #can also change to path (e.g. 'Sub Directory/graph.png') to set location for saved plot)