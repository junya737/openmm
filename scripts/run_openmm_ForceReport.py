from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


class ForceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        for f in forces:
            self._out.write('[%g %g %g],' % (f[0], f[1], f[2]))


# Create PDBFile object
pdb = PDBFile('../examples/input.pdb')
# ForceField object
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# System object
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)


# Integrator object that propagates the dynamics
# thetemperature (300 K), the friction coefficient (1 ps-1), and the step size (0.004 ps)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
# combines the molecular topology, system, and integrator to begin a new simulation.
simulation = Simulation(pdb.topology, system, integrator)

# Specifies the initial atom positions
simulation.context.setPositions(pdb.positions)
# perform a local energy minimization
simulation.minimizeEnergy()
# creates a “reporter” to generate output, and the timestep is 1000
simulation.reporters.append(PDBReporter('output.pdb', 1))
simulation.reporters.append(ForceReporter('forces.txt', 1))

# adds another reporter to print out some basic information every 1000 time steps, and stdout means the output will be printed to the screen.
simulation.reporters.append(StateDataReporter(stdout, 1, step=True, potentialEnergy=True, temperature=True))
simulation.step(3)