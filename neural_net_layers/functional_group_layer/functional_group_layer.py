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

common_functional_groups_smarts = {
    'imidazole': '[#6]1:[#6]:[#7]:[#6]:[#7H]:1',
    'halogen': '[F,Cl,Br,I]',
    'acetic anhydride':'[CX3](=[OX1])[OX2][CX3](=[OX1])',
    'hydroxyl': '[OX2H]', 'phenol': '[OX2H][cX3]:[c]',
    'aldehyde': '[CX3H1](=O)[#6]',
    'amide': '[NX3][CX3](=[OX1])[#6]',
    'ketone': '[#6][CX3](=O)[#6]',
    'ester': '[#6]-C(=O)O-[#6]',
    'ether': '[OD2]([#6])[#6]',
    'amine': '[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]',
    'enamine': '[NX3][CX3]=[CX3]',
    'vinylic alkene': '[$([CX3]=[CX3])]',
    'alkyne': '[$([CX2]#C)]',
    'allenic alkene': '[$([CX2](=C)=C)]',
    'sulfide': '[#16X2H0]',
    'sulfonamide': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]',
    'sulfoxide': '[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]',
    'thiazole': '[#6]1:[#6]:[#16]:[#6]:[#7]:1',
    'cyclopropane': '[#6]1-[#6]-[#6]-1'
}

rings_in_drugs_smarts = {
    'benzene': '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1',
    'pyridine': '[#6]1:[#6]:[#6]:[#6]:[#7]:[#6]:1',
    'piperidine': '[#7]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
    'piperazine': '[#7]1-[#6]-[#6]-[#7]-[#6]-[#6]-1',
    'cyclohexane': '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
    'oxane': '[#8]1-[#6]-[#6]-[#6]-[#6]-[#6]-1',
    'imidazole': '[#6]1:[#7]:[#6]:[#6]:[#7H]:1',
    'pyrrolidine': '[#6]1-[#6]-[#6]-[#7]-[#6]-1',
    '(R)-5-thia-1-azabicyclo[4.2.0]oct-2-en-8-one': '[#8]=[#6]1-[#6]-[#6@@H]2-[#7]-1-[#6]=[#6]-[#6]-[#16]-2',
    'cyclopropane': '[#6]1-[#6]-[#6]-1',
    'tetrahydrofuran': '[#6]1-[#6]-[#6]-[#8]-[#6]-1',
    'thiazole': '[#6]1:[#7]:[#6]:[#6]:[#16]:1',
    'indole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#7H]:2',
    'diazine': '[#6]1:[#7]:[#6]:[#6]:[#6]:[#7]:1',
    '(R)-4-thia-1-azabicyclo[3.2.0]heptan-7-one': '[#8]=[#6]1-[#7]2-[#6]-[#6]-[#16]-[#6@@H]-2-[#6]-1',
    '4-((3aS)-octahydro-1H-inden-5-yl)-3-propylcyclohexa-2,5-dien-1-one': '[#8]=[#6]1-[#6]=[#6]-[#6](-[#6]2-[#6]-[#6@@H]3-[#6]-[#6]-[#6]-[#6]-3-[#6]-[#6]-2)-[#6](-[#6]-[#6]-[#6])=[#6]-1',
    'tetrazole': '[#7]1=[#7]-[#7]=[#6]-[#7]-1',
    'cyclopentane': '[#6]1-[#6]-[#6]-[#6]-[#6]-1',
    'thiophenyl': '[#6]1:[#6]:[#6]:[#6]:[#16]:1',
    'naphthalene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
    '1H-benzo[d]imidazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#7]:[#6]:[#7H]:2',
    'quinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:[#7]:2',
    '1H-purine': '[#6]12:[#6]:[#7H]:[#6]:[#7]:[#6]-1:[#7]:[#6]:[#7]:2',
    '4-((3aS,5S)-octahydro-1H-inden-5-yl)-3-propylcyclohex-2-en-1-one': '[#8]=[#6]1-[#6]-[#6]-[#6](-[#6@@H]2-[#6]-[#6@@H]3-[#6]-[#6]-[#6]-[#6]-3-[#6]-[#6]-2)-[#6](-[#6]-[#6]-[#6])=[#6]-1',
    'furan': '[#6]1:[#6]:[#6]:[#6]:[#8]:1',
    '1H-1,2,4-Triazole': '[#7]1:[#6]:[#7]:[#6]:[#7H]:1',
    '10H-Phenothiazine': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7]-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#16]-2',
    'quinazoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#6]:[#7]:2',
    'morpholine': '[#6]1-[#6]-[#7]-[#6]-[#6]-[#8]-1',
    'pyrimidin-2(1H)-one': '[#8]=[#6]1:[#7]:[#6]:[#6]:[#6]:[#7H]:1',
    'quinolin-4(1H)-one': '[#8]=[#6]1:[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2):[#7H]:[#6]:[#6]:1',
    '4-((3aS,5S)-octahydro-1H-inden-5-yl)-3-propylcyclohexa-2,5-dien-1-one': '[#8]=[#6]1-[#6]=[#6]-[#6](-[#6@H]2-[#6]-[#6]-[#6]3-[#6@@H](-[#6]-[#6]-[#6]-3)-[#6]-2)-[#6](-[#6]-[#6]-[#6])=[#6]-1',
    'isoxazole': '[#6]1:[#6]:[#6]:[#7]:[#8]:1',
    'imidazoline': '[#6]1=[#7]-[#6]-[#6]-[#7]-1',
    '1,4-dihydropyridin': '[#6]1=[#6]-[#6]-[#6]=[#6]-[#7]-1',
    'pyrimidine-2,4(1H,3H)-dione': '[#8]=[#6]1:[#7H]:[#6](:[#6]:[#6]:[#7H]:1)=[#8]',
    '3,4-dihydro-2H-benzo[e][1,4]diazepin-2-one': '[#8]=[#6]1-[#7]=[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2=[#6]-[#7]-[#6]-1',
    'cyclohexene': '[#6]1=[#6]-[#6]-[#6]-[#6]-[#6]-1',
    'pyrrolidin-2-one': '[#8]=[#6]1-[#7]-[#6]-[#6]-[#6]-1',
    'imidazolidine-2,4-dione': '[#8]=[#6]1-[#6]-[#7]-[#6](-[#7]-1)=[#8]',
    '1,2,3,4-tetrahydroisoquinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#7]-[#6]-2',
    '3,4-dihydro-2H-benzo[e][1,2,4]thiadiazine 1,1-dioxide': '[#8]=[#16]1(-[#7]-[#6]-[#7]-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#8]',
    '(3aS,5S)-5-(2-propylphenyl)octahydro-1H-indene': '[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6@@H]1-[#6]-[#6@@H]2-[#6]-[#6]-[#6]-[#6]-2-[#6]-[#6]-1',
    '1H-pyrazole': '[#7]1:[#6]:[#6]:[#6]:[#7H]:1',
    'quinuclidine': '[#6]12-[#6]-[#6]-[#7](-[#6]-[#6]-1)-[#6]-[#6]-2',
    'epoxide': '[#6]1-[#6]-[#8]-1',
    'pyrazine': '[#6]1:[#6]:[#7]:[#6]:[#6]:[#7]:1',
    'oxazolidinone': '[#8]=[#6]1-[#8]-[#6]-[#6]-[#7]-1',
    'tetrahydronaphthalene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]-[#6]-2',
    'adamantane': '[#6]12-[#6]-[#6]3-[#6]-[#6](-[#6]-1)-[#6]-[#6](-[#6]-3)-[#6]-2',
    '1,8-naphthyridin-4(1H)-one': '[#8]=[#6]1:[#6]:[#6]:[#7H]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#7]:2',
    '3,7-dihydro-1H-purine-2,6-dione': '[#8]=[#6]1:[#6]2:[#7H]:[#6]:[#7]:[#6]:2:[#7H]:[#6](:[#7H]:1)=[#8]',
    '(3aS,5S)-5-((2S)-2-propylcyclohexyl)octahydro-1H-indene': '[#6]-[#6]-[#6]-[#6@H]1-[#6]-[#6]-[#6]-[#6]-[#6]-1-[#6@@H]1-[#6]-[#6@@H]2-[#6]-[#6]-[#6]-[#6]-2-[#6]-[#6]-1',
    '7,8,9,10-tetrahydrotetracene-5,12-dione': '[#8]=[#6]1-[#6]2:[#6]:[#6]3-[#6]-[#6]-[#6]-[#6]-[#6]:3:[#6]:[#6]:2-[#6](=[#8])-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2',
    'cyclobutane': '[#6]1-[#6]-[#6]-[#6]-1',
    '1,2-dihydro-3H-1,2,4-triazol-3-one': '[#8]=[#6]1:[#7H]:[#7H]:[#6]:[#7]:1',
    '1,3,4-thiadiazole': '[#6]1:[#7]:[#7]:[#6]:[#16]:1',
    'azepane': '[#6]1-[#7]-[#6]-[#6]-[#6]-[#6]-[#6]-1',
    '8-azabicyclo[3.2.1]octane': '[#6]12-[#6]-[#6]-[#6]-[#6](-[#6]-[#6]-1)-[#7]-2',
    'piperidine-2,6-dione': '[#8]=[#6]1-[#7]-[#6](-[#6]-[#6]-[#6]-1)=[#8]',
    '2,3-dihydro-1H-indene': '[#8]=[#6]1-[#7]-[#6](-[#6]-[#6]-[#6]-1)=[#8]',
    'benzo[d]isoxazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#8]:2',
    '1,9-dihydro-6H-purin-6-one': '[#8]=[#6]1:[#6]2:[#6](:[#7H]:[#6]:[#7]:2):[#7]:[#6]:[#7H]:1',
    '9H-fluorene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#6]-2',
    '10,11-dihydro-5H-dibenzo[b,f]azepine': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#7]-2',
    '(6aR,10aR)-4,6,6a,7,8,9,10,10a-octahydroindolo[4,3-fg]quinoline': '[#6]12:[#6]:[#6]:[#6]:[#6]3:[#6]:1:[#6](-[#6]-[#6@@H]1-[#6@@H]-2-[#6]-[#6]-[#6]-[#7]-1):[#6]:[#7H]:3',
    '1H-pyrrole': '[#6]1:[#6]:[#6]:[#6]:[#7H]:1',
    '1,3-dioxolane': '[#8]1-[#6]-[#6]-[#8]-[#6]-1',
    '(1R,5S)-3-azabicyclo[3.1.0]hexane': '[#6]1-[#6@@H]2-[#6@H]-1-[#6]-[#7]-[#6]-2',
    'cyclopentanone': '[#8]=[#6]1-[#6]-[#6]-[#6]-[#6]-1',
    'pyrrolidine-2,5-dione': '[#8]=[#6]1-[#7]-[#6](-[#6]-[#6]-1)=[#8]',
    'pyrazolidine': '[#8]=[#6]1-[#7]-[#7]-[#6](-[#6]-1)=[#8]',
    '(R)-1-azabicyclo[3.2.0]hept-2-en-7-one': '[#8]=[#6]1-[#7]2-[#6]=[#6]-[#6]-[#6@@H]-2-[#6]-1',
    'thiazolidine-2,4-dione': '[#8]=[#6]1-[#6]-[#16]-[#6](-[#7]-1)=[#8]',
    'benzofuran': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#8]:2',
    '1H-indazole': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#7]:[#7H]:2',
    'indolin-2-one': '[#8]=[#6]1-[#7]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6]-1',
    'benzo[b]thiophene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#16]:2',
    '(R)-1,2,3,7,8,8a-hexahydronaphthalene': '[#6]12=[#6]-[#6]-[#6]-[#6]-[#6@@H]-1-[#6]-[#6]-[#6]=[#6]-2',
    '4,5,6,7-tetrahydrothieno[3,2-c]pyridine': '[#6]12:[#6]:[#6]:[#16]:[#6]:1-[#6]-[#6]-[#7]-[#6]-2',
    '4H-chromen-4-one': '[#8]=[#6]1:[#6]:[#6]:[#8]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2',
    '3,4-dihydroquino-2(1H)-one': '[#8]=[#6]1-[#6]-[#6]-[#6]2:[#6](-[#7]-1):[#6]:[#6]:[#6]:[#6]:2',
    'napthalene-1,4-dione': '[#8]=[#6]1-[#6]=[#6]-[#6](=[#8])-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2',
    '2H-benzo[e][1,2,4]thiadiazine 1,1-dioxide': '[#8]=[#16]1(-[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2)-[#7]=[#6]-[#7]-1)=[#8]',
    '4H-benzo[f][1,2,4]triazolo[4,3-a][1,4]diazepine': '[#6]12-[#7]3:[#6](-[#6]-[#7]=[#6]-[#6]:1:[#6]:[#6]:[#6]:[#6]:2):[#7]:[#7]:[#6]:3',
    '9H-thioxanthene': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]1:[#6](:[#6]:[#6]:[#6]:[#6]:1)-[#16]-2',
    'dummy_37': '[#8]=[#6]1-[#8]-[#6]-[#6@H]2-[#6@H]-1-[#6]-[#6]1:[#6](-[#6]-2):[#6]:[#6]2:[#6](:[#6]:1)-[#8]-[#6]-[#8]-2',
    'dummy_39': '[#6]12:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6@@]13-[#6@H](-[#6]-[#6]-[#6]4-[#6@@H]-1-[#7](-[#6]-[#6]=[#6]-4)-[#6]-[#6]-3)-[#7]-2',
    'dummy_40': '[#8]=[#6]1-[#6]2=[#6]-[#6]3-[#6@H](-[#6]-[#6]=[#6]-[#6]-3=[#8])-[#6]-[#6@@H]-2-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1',
    '1H-1,2,3-triazole': '[#7]1:[#7]:[#6]:[#6]:[#7H]:1',
    'azetidin-2-one': '[#8]=[#6]1-[#7]-[#6]-[#6]-1',
    'oxetan-2-one': '[#8]=[#6]1-[#8]-[#6]-[#6]-1'
}

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

    for key, value in common_functional_groups_smarts.items():

        try:
            if Chem.MolFromSmarts(value) is not None:
                validated_smarts.append([key, value, Chem.MolFromSmarts(value)])
        except:
            pass

    for key, value in rings_in_drugs_smarts.items():

        try:
            if Chem.MolFromSmarts(value) is not None:
                validated_smarts.append([key, value, Chem.MolFromSmarts(value)])
        except:
            pass



    for compound in compounds:

        smiles = compound['smiles']
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
        molecule = Chem.MolFromSmiles(smiles)

        matches = []

        for i in range(0, len(validated_smarts)):

            row = validated_smarts[i]
            pattern_name = row[0]
            pattern_mol = row[2]

            try:
                substructs = molecule.GetSubstructMatches(pattern_mol)
                if substructs:
                    matches.append(pattern_name)

            except:
                pass

        compound['matches'] = matches

    submatches = {}

    for compound in compounds:

        matches = compound['matches']

        for match in matches:
            relations = [match + ' - ' + i for i in matches if i != match]

            for relation in relations:
                if relation not in submatches:
                    submatches[relation] = 1
                else:
                    submatches[relation] += 1

            compound['relations'] = relations

    functional_group_network_data = {}

    for compound in compounds:

        for i in compound['matches']:

            if i not in functional_group_network_data:

                functional_group_network_data[i] = {}


    for compound in compounds:

        if 'relations' in compound:

            for i in compound['relations']:

                match = i.split(' - ')[0]
                relation = i.split(' - ')[1]

                if relation in functional_group_network_data[match]:

                    functional_group_network_data[match][relation] += 1

                else:

                    functional_group_network_data[match][relation] = 0

    visualization_network = []

    for key, value in functional_group_network_data.items():

        row = {}

        row['name'] = key
        row['size'] = list(value.values())
        row['imports'] = list(value.keys())

        visualization_network.append(row)

    with open('functionaL_group.json', 'w') as outfile:
         json.dump(visualization_network, outfile)
