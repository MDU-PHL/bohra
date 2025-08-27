#!/usr/bin/env python3

import sys, json, pathlib


available_species = json.load(open(f"{pathlib.Path(__file__).parent.resolve() / 'species.json'}", 'r'))
species_obs = sys.argv[1].strip()

if species_obs.replace(" ", "_") in available_species['abritamr_species'] :
    print(f"-sp {species_obs.replace(' ', '_')}")
elif species_obs.split()[0] in available_species['abritamr_species'] :
    print(f"-sp {species_obs.split()[0]}")
elif species_obs in available_species['abritamr_synonyms'] :
    print(f"-sp {available_species['abritamr_synonyms'][species_obs]}")
else:
    print("")