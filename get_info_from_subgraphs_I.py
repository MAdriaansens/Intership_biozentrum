#with this script we used a known dataset so we had some candidates for best protein hits etc, so these were manually imported.
#as imput we used data obtained from graphs of networkx

import os
import networkx
import sys
import subprocess
import argparse
import glob
import json
from subprocess import check_output
import re

cluster_percentage = sys.argv[1]


# orginal path results/Clustered_100/subgraphs/
hits = list(glob.glob('results/Clustered_{}/subgraphs/*.gml'.format(cluster_percentage)))

#get UniprotDB
sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')



MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

from src import extract_uniprot  as uniprot
uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

from src import extract_interpro  as interpro
interpro_db = interpro.interpro_db_diggested(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

#set up empty folders

output_list_prot = []
output_list_sg = []
protein_ids_dict = {}


def is_it_mixed(hits, cluster_percentage):

    count = 0
    output_uniprot =[]
    temp_protein_Id = []
    outfolder = './results/Clustered_{}/json_uniprot_id'.format(cluster_percentage)

    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    for hit in hits:
        subgraph_id = hit.split("_")[-1].split(".")[0]
        text = subprocess.check_output(['grep', 'label', '-h', '{}'.format(hit)], universal_newlines = True)

        protein_id_devided = re.findall(pattern,text)
        uniprot_id = protein_id_devided[0]
        output_uniprot.append(uniprot_id)

    
        with open('{}/uniprot_ids_{}.json'.format(outfolder, subgraph_id), 'w') as f:
            json.dump(output_uniprot, f)



def search_uniprot(uniprot_db, output_list_prot):
    subgraphs = list(glob.glob('results/Clustered_{}/json_uniprot_id/*.json'.format(cluster_percentage)))
    for subgraph in subgraphs:
        subgraph_id = subgraph.split(".")[0]
        
        with open('{}'.format[subgraph]) as json_file:
            


    list_protein_function = []

    list_protein_kingdom = []

    target_ACCs = output_list_prot #list of uniprot ids
    curr_collection= uniprot_db.col.find({ '_id': { "$in": target_ACCs }})
    Bacteria_count = 0
    Archaea_count = 0
    Eukaryota_count = 0
    Else_count = 0

    for document in curr_collection:
        #curr_name is dictionairy
        #curr_taxid is list of lists

        curr_acc = document['_id']
        curr_data = document['data']
        curr_taxid = curr_data['TAXID']
        curr_name = curr_data['NAME']


        try:
            key_1 =  list(curr_name)[0]
            protein_name = curr_name.get(key_1)
            list_protein_function.append(protein_name)
        except IndexError:
            print('not enough keys')



        flat_tax = [v for item in curr_taxid for v in (item if isinstance(item,list) else [item])]
        if 'Bacteria' in flat_tax:
            kingdom = 'Bacteria'
            Bacteria_count = Bacteria_count + 1
        elif 'Eukaryota' in flat_tax:
            kingdom = 'Eukaryota'
            Eukaryota_count = Eukaryota_count +1
        elif 'Archaea' in flat_tax:
            kingdom = 'Archaea'
            Archaea_count = Archaea_count + 1
        else:
            kingdom = 'non Archaea,Bacteria or Eukarya'
            Else_count = Else_count + 1

        #list_protein_kingdom = list_protein_kingdom.append(kingdom)
#we became very interested in what our proteins contained so via extraction of a list of the most common function we performed a count of each protein. 
#This would later be removed in a second iteration of this code using the max() functon + count
        count_KH = (sum('KH' in s for s in list_protein_function))
        count_NusA = (sum('NusA' in s for s in list_protein_function))
        count_Ribonuclease = (sum('Ribonuclease' in s for s in list_protein_function))
        count_Flagellar = (sum('Flagellar' in s for s in list_protein_function))
        count_ATPase = (sum('ATPase' in s for s in list_protein_function))
        count_30S = (sum('30S' in s for s in list_protein_function))
        count_KhpA = (sum('KhpA' in s for s in list_protein_function))
        count_GTPase = (sum('GTPase' in s for s in list_protein_function))
        count_R3H = (sum('R3H' in s for s in list_protein_function))
        count_SpoIIIJ = (sum('SpoIIIJ' in s for s in list_protein_function))
        count_Uncharacterized = (sum('Uncharacterized protein' in s for s in list_protein_function))
        count_PhoH = (sum('PhoH' in s for s in list_protein_function))
        count_Betalactamase = (sum('Beta-lactamase' in s for s in list_protein_function))
        count_Branchpoint = (sum('Branchpoint' in s for s in list_protein_function))
        count_S1_motif = (sum('S1 motif' in s for s in list_protein_function))


   # print('\n' + 'Distribution of kingdoms:' + '\n' + 'no of Bacteria: {}'.format(Bacteria_count) + '\n' + 'no of Eukaryota: {}'.format(Eukaryota_count) + '\n'+ 'no of Archaea: {}'.format(Archaea_count) + '\n' + 'no of Else: {}'.format(Else_count) + '\n')

  #  print('\n' + 'Distribution of functions:' + '\n' + 'No KH_domain: {}'. format(count_KH) + '\n' + 'No S1_motif: {}'.format(count_S1_motif) + '\n' + 'No ATPase: {}'.format(count_ATPase) \
   # + '\n' + 'No KhpA: {}'.format(count_KhpA) + '\n' + 'No Beta_lactamase: {}'.format(count_Betalactamase) + '\n' + 'No GTPase: {}'.format(count_GTPase) + '\n' + 'No Ribonuclease: {}'.format(count_Ribonuclease) +'\n' + 'No Uncharacterized proteins: {}'.format(count_Uncharacterized) + '\n' + 'No SpoIIIJ: {}'.format(count_SpoIIIJ) + '\n' + 'No: 30S {}'.format(count_30S) + '\n' + 'No Flagellar: {}'.format(count_Flagellar)\
   # + '\n' + 'No PhoH: {}'.format(count_PhoH) + '\n' + 'No R3H: {}'.format(count_R3H)  + '\n' + 'No NusA: {}'.format(count_NusA))


#maincode
search_uniprot(uniprot_db, output_list_prot)
is_it_mixed(hits, cluster_percentage)
                                       
