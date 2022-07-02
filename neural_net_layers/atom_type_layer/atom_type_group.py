# imports
# -------
import os
import json

# tinydbs
# -------
from tinydb import TinyDB, Query

# rdkit
# -----
from rdkit import Chem

user = Query()

if __name__ == '__main__':

    database_path = '../../tinydbs'
    _, _, tinydb_paths = next(os.walk('../../tinydbs/'))

    tinydbs = []

    for i in tinydb_paths[0:8]:

        db = TinyDB(os.path.join(database_path, i))
        compounds = db.search(user.smiles.exists())

        tinydbs.append(compounds)

        print ("Loaded: %s" % i)

    compounds = sum(tinydbs, [])
    validated_smarts = []


    for compound in compounds:

        if 'atom_penalties' in compound:
            atoms = compound['atom_penalties']
        else:
            atoms = compound['charge_atom_penalties']

        atom_types = []

        for i in atoms:
            atom_types.append(i[0])

        compound['atom_types'] = atom_types


    submatches = {}

    for compound in compounds:

        matches = compound['atom_types']

        for match in matches:
            relations = [match + ' - ' + i for i in matches if i != match]

            for relation in relations:
                if relation not in submatches:
                    submatches[relation] = 1
                else:
                    submatches[relation] += 1

            compound['atom_type_relations'] = relations

    atom_type_network_data = {}

    for compound in compounds:

        for i in compound['atom_types']:

            if i not in atom_type_network_data:

                atom_type_network_data[i] = {}


    for compound in compounds:

        if 'atom_type_relations' in compound:

            for i in compound['atom_type_relations']:

                match = i.split(' - ')[0]
                relation = i.split(' - ')[1]

                if relation in atom_type_network_data[match]:

                    atom_type_network_data[match][relation] += 1

                else:

                    atom_type_network_data[match][relation] = 0

    visualization_network = []

    for key, value in atom_type_network_data.items():

        row = {}

        row['name'] = key
        row['size'] = list(value.values())
        row['imports'] = list(value.keys())

        visualization_network.append(row)

    with open('atom_type_group.json', 'w') as outfile:
        json.dump(visualization_network, outfile)