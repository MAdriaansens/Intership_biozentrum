#way to get information from protein_ids of the seperate subgraphs using searches against Interpro. Also included kingdom origin however is excluded due to focus mainly being on Archaea specifc proteins
#follow up step is integration of multiprocessing to speed it up

import os
import networkx as nx
import sys
import argparse
import glob
import json
import re


infolder = sys.argv[1]
outfolder = '{}/json_uniprot_id'.format(infolder)
print(infolder)
print(outfolder)

hits = list(glob.glob('{}/subgraphs/*.gml'.format(infolder)))


subgraphs = list(glob.glob('{}/subgraphs/*json'.format(infolder)))

#get UniprotDB
sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')
from src import extract_interpro  as interpro
from src import extract_uniprot  as uniprot

MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
interpro_db = interpro.interpro_db_diggested(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

def is_it_mixed(hits,outfolder):
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    for hit in hits:
        subgraph_id = hit.split("_")[-1].split(".")[0]
        file1 = open(hit, 'r')
        pattern = r'(?<=\")(.*?)(?=\{)'

        new_list = []
        for line in file1:
            parsed = re.findall(pattern, line)
            ##print(parsed)
            if parsed != []:
                temp_protein_id = parsed[0]
                new_list.append(temp_protein_id)
                #print("these protein ids belong to subgraph: {}".format(subgraph_id), temp_protein_id)

                with open('{}/protein_id_subgraph_{}.json'.format(outfolder, subgraph_id), 'w') as f:
                   json.dump(new_list, f)
subgraphs = list(glob.glob('{}/*.json'.format(outfolder)))

def search_kingdom_vs_uniprot(uniprot_db, subgraphs, outfolder):
    #counts the kingdoms but also count the function of the most common protein
    #count_eukarya_only = 0
    #count_archaea_only = 0
    #count_prokarya_mixed = 0
    #count_bacteria_only = 0
    #count_mixed = 0
    lists_values = []
    count = 0
    if not os.path.isdir("{}/protein_information".format(outfolder)):
        os.mkdir("{}/protein_information".format(outfolder))

    for subgraph in subgraphs:
        subgraph_id = subgraph.split("_")[-1].split(".")[0]
        protein_information = []

        #count if needed
        mixed_count = 0
        Bacteria_count = 0
        Archaea_count = 0
        Eukaryota_count = 0
        Else_count = 0

        target_ACCs = json.load(open(subgraph,'r'))
        for i in target_ACCs:
            curr_collection= interpro_db.query('{}'.format(i))
            for document in curr_collection:
                curr_acc = document['_id']
                curr_data = document['data']
                curr_name = curr_data[-1]
                curr_name = ' '.join(str(e) for e in curr_name)
                curr_name = curr_name.split(",")[0].split("[")[0]

                #subgraph_info = ['kingdom distribution for subgraph {}:'.format(subgraph_id) + "Bacteria:{} ".format(Bacteria_count) + 'Eukarya:{} '.format(Eukaryota_count) +'Archaea:{} '.format(Archaea_count)+ 'Else:{} '.format(Else_count)]

                inform = [curr_name] #i and kingdom can be added to give these values as well.
                protein_information.append(inform)

        mostFrequent = max(protein_information, key=protein_information.count, default ="no_protein_id")
        subgraph_info_json = []
        subgraph_info_json.append(protein_information)
        with open('{}/protein_information/information_protein_id_subgraph_{}.json'.format(outfolder, subgraph_id), 'w') as f:
            json.dump(subgraph_info_json, f)

        if mostFrequent == "no_protein_id":
            count = count + 1
    print(count)

is_it_mixed(hits,outfolder)
search_kingdom_vs_uniprot(uniprot_db, subgraphs, outfolder)

