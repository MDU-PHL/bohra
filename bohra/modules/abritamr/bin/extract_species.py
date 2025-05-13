#!/usr/bin/env python3

import sys, json


available_species = json.load(open(sys.argv[1], 'r'))
species_obs = sys.argv[2]

if species_obs.replace(" ", "_") in available_species['abritamr_species'] :
    print(f"-sp {species_obs.replace(' ', '_')}")
elif species_obs.split()[0] in available_species['abritamr_species'] :
    print(f"-sp {species_obs.split()[0]}")

else:
    print("")