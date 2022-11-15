##################################################################################
'''
Step zero of the vibronic model workflow.

Assigns constants and values used across multiple steps. Records this info
currently in 'master.py' file; anticipated to be replaced by JSON file in the future.

'''
###################################################################################

# system imports
import json

# ask the user for how much memory and processors desired
molecule_name = input("Name of molecule: ")
natoms = input("Amount of atoms in molecule: ")
nstates = input("Amount of excited states to calculate: ")
ndisplacements = input("Desired number of displacements along each normal mode, e.g. 2: ")
nfreqs = int(3*int(natoms)-6)

# store data in master file
with open('master_values.json', 'w') as fp:

    dictionary = {
        "molecule_name": molecule_name,
        "num_atoms": int(natoms),
        "num_excited_states": int(nstates),
        "num_freqs": int(nfreqs),
        "num_displacements": int(ndisplacements),
        "scaling_factor": '1.d0',
        "sorting_atoms": 'no',
        "step1": {
            "nproc":8,
            "mem": 30
            },
        "step1b": {
            "nproc": 7,
            "mem": 25
            },
        "step4": {
            "nproc": 10,
            "mem": 35
            },
        "step5": {
            "nproc": 11,
            "mem": 40
            }
    }

    json.dump(dictionary, fp, indent=0)
