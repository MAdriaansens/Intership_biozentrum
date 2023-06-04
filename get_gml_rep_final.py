import os
import networkx as nx
import pickle
import sys
import subprocess
import argparse
import glob
import json
import numpy as np
import re

cluster_percentage = sys.argv[1]
outfolder = './results/Clustered_{}/json_uniprot_id'.format(cluster_percentage)

# orginal path results/Clustered_100/subgraphs/
hits = list(glob.glob('results/Clustered_{}/subgraphs/*.gml'.format(cluster_percentage)))
subgraphs = list(glob.glob('{}/*json'.format(outfolder)))

#get UniprotDB
sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')
from src import extract_interpro  as interpro 
from src import extract_uniprot  as uniprot
 
MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

 
uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
interpro_db = interpro.interpro_db_diggested(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

def parse_gml(hits,outfolder):
      
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    for hit in hits:
        subgraph_id = hit.split("_")[-1].split(".")[0]
        file1 = open(hit, 'r')
        pattern = r'(?<=\_)(.*?)(?=\_)'
        #print(pattern)
        new_list = []
        for line in file1:
            parsed = re.findall(pattern, line)
            if parsed != []:
                temp_protein_id = parsed[0]
                new_list.append(temp_protein_id)
                #print("these protein ids belong to subgraph: {}".format(subgraph_id), temp_protein_id)
                with open('{}/protein_id_subgraph_{}.json'.format(outfolder, subgraph_id), 'w') as f:
                    json.dump(new_list, f)

def find_kingdom_name(uniprot_db, interpro_db, subgraphs, outfolder):
    output_protein_info = '{}/protein_information'.format(outfolder)
    if not os.path.isdir(output_protein_info):
        os.mkdir(output_protein_info)
    count_eukarya_only = 0
    count_archaea_only = 0
    count_prokarya_mixed = 0
    count_bacteria_only = 0
    count_mixed = 0    
    lists_values = []
    subgraphs = list(glob.glob('{}/*json'.format(outfolder)))
    for subgraph in subgraphs:      
        

        
        subgraph_info_json = []
        protein_information = []
        
        mixed_count = 0
        Bacteria_count = 0
        Archaea_count = 0
        Eukaryota_count = 0
        Else_count = 0

        subgraph_id = subgraph.split("_")[-1].split(".")[0]
        target_ACCs = json.load(open(subgraph,'r'))
        for i in target_ACCs:
 #        curr_collection= uniprot_db.col.find({ '_id': { "$in": target_ACCs }})[0]
            curr_collection= uniprot_db.query('{}'.format(i))
            #print(curr_collection[0])

            for document in curr_collection:
           #    curr_name is dictionairy
           #    curr_taxid is list of lists
               
               curr_acc = document['_id']
               curr_data = document['data']
               curr_taxid = curr_data['TAXID']
               curr_name = curr_data['NAME']
               
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
               
               curr_collection= interpro_db.query('{}'.format(i))
               for document in curr_collection:
                   curr_acc = document['_id']
                   curr_data = document['data'] 
                   curr_name = curr_data[-1] 
                   curr_name = ' '.join(str(e) for e in curr_name)
                   curr_name = curr_name.split(",")[0].split("[")[0] 
                   subgraph_info = ['kingdom distribution for subgraph {}:'.format(subgraph_id) + "Bacteria:{} ".format(Bacteria_count) + 'Eukarya:{} '.format(Eukaryota_count) +'Archaea:{} '.format(Archaea_count)+ 'Else:{} '.format(Else_count)]
                   inform = [curr_name] #i and kingdom can be added to give these values as well. 
                   protein_information.append(inform)
        
         
        mostFrequent = max(protein_information, key=protein_information.count)
        print(subgraph_id, mostFrequent) 
        
        if Else_count != 0:
            print("else!!")
            break
        else:

            if Bacteria_count != 0 and Archaea_count != 0 and Eukaryota_count != 0:
                #print(subgraph_id + ' is mixed')
                count_mixed = count_mixed + 1       
            if Bacteria_count != 0 and Archaea_count == 0 and Eukaryota_count == 0:
                #print(subgraph_id + ' is bacteria only')
                count_bacteria_only = count_bacteria_only + 1
            if Bacteria_count == 0 and Archaea_count == 0 and Eukaryota_count != 0:
                #print(subgraph_id + ' is eukarya only')
                count_eukarya_only = count_eukarya_only + 1
            if Bacteria_count == 0 and Archaea_count == 0 and Eukaryota_count != 0:
                #print(subgraph_id + ' prokarya-mixed')
                count_prokarya_mixed = count_prokarya_mixed + 1
            if Bacteria_count == 0 and Archaea_count != 0 and Eukaryota_count == 0:
               # print(subgraph_id + ' archaea only')
                count_archaea_only = count_archaea_only + 1
        
                    
         #put all subgraph info in seperate files per subgraph           
        subgraph_info_json.append(subgraph_info)
        subgraph_info_json.append(protein_information)
 #       print(subgraph_info_json)
        with open('{}/information_protein_id_subgraph_{}.json'.format(output_protein_info, subgraph_id), 'w') as f:
            json.dump(subgraph_info_json, f)
           
    #print("at a {} percentage coverage the kingdoms consistency of the connected components consists of the following:".format(cluster_percentage) + '\nBacteria only {}'.format(count_bacteria_only) + '\nEukaryota only {}'.format(count_eukarya_only) \
            #         + '\nmixed {}'.format(count_mixed) + '\nArchaea only {}'.format(count_archaea_only) + '\nProkaryota only {}'.format(count_prokarya_mixed))   

parse_gml(hits,outfolder)
find_kingdom_name(uniprot_db, interpro_db, subgraphs, outfolder)
