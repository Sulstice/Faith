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

    db = '/Users/sulimansharif/mackerell_group/faith/global_chem/atom_type_group_new.json'
    db = TinyDB(db)

    validated_smarts = []
    compounds = db.search(user.smiles.exists())

    for compound in compounds:

        atom_types = [row.split()[2] for row in compound['atom_rows'] ]
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

    with open('atom_type_group_new.json', 'w') as outfile:
        json.dump(visualization_network, outfile)