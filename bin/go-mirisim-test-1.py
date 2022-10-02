#!/usr/bin/env python

from mirisim import MiriSimulation
from mirisim.config_parser import SimConfig, SceneConfig, SimulatorConfig

#MiriSimulation.generate_configfiles()

#mysim = MiriSimulation.from_configfiles('simulation.ini')

simcfg = SimConfig('mrs_simulation.ini')
scenecfg = SceneConfig('scene.ini')
simulatorcfg = SimulatorConfig('simulator.ini')
mysim = MiriSimulation(simcfg, scenecfg, simulatorcfg)

mysim.run()
