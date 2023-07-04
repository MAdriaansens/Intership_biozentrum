#open subgraph gml, make option for a list as well as this is the representatives
import glob
import os
import math
import sys
import networkx as nx
from networkx.algorithms.centrality import closeness_centrality


subgraphs = list(glob.glob('results/test_gml/*gml'))

sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')
from src import extract_interpro  as interpro
from src import extract_uniprot  as uniprot

MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

def get_centrality(subgraphs):
    
    #goal of this function is to find the nodes with the most clossest centrality of all the nodes, it skipps subgraphs below 3 nodes it also looks if there are more than one clossest node but still picks a random one. (maybe include all of these anyhow and then pick longest one
    subgraph_id_closest_node = {}
    list_longer_subgraph = []
    length_dict = {}
    top_percent_cent = {}
    dict_representatives = {}
    
    for subgraph in subgraphs:
        count = 0
        rep_potential = []
        subgraph_id = subgraph.split("_")[-1].split(".")[0]
        G = nx.read_gml(subgraph)
        closeness_subgraph = closeness_centrality(G)
        #this is to skipp subgraphs containing only 1 or 2 nodes as centrality or representatives between these are not needed or that different
        if len(closeness_subgraph) <3:
            pass
        else:
            length_dict = get_length(closeness_subgraph, uniprot)
            if length_dict:

                #print(length_dict)
                #list subgraphs with more than 2 values
                list_longer_subgraph.append(subgraph_id)
                highest_closeness = max(closeness_subgraph, key=closeness_subgraph.get)


                #we also want to make a list of the top most central proteins, right now 25% is taken as the number of nodes in chosen subgraphs is between 4 and 12
                N = math.ceil((len(closeness_subgraph) * 0.25))
                top_percent_cent = dict(sorted(closeness_subgraph.items(), key = lambda x: x[1], reverse = True)[:N])
                
                #only take values exceeding the average 
                res = 0
                for val in closeness_subgraph.values():
                        res += val
                avg = res/len(closeness_subgraph)
                
                for key, value  in top_percent_cent.items():
                    centrality_devided_byAvg = value/avg
                    top_percent_cent[key] = centrality_devided_byAvg
                   # print(top_percent_cent)
                        #print('centrality', key, centrality_devided_byAvg)
            else:
                print('empty')      

            rep_pres = potential_representative(length_dict, top_percent_cent)
            weighted_values = multiply(top_percent_cent, length_dict)
            if weighted_values:
                representative = (max(weighted_values, key=weighted_values.get))
                dict_representatives[subgraph_id] = representative
    print(dict_representatives)
            

            
#now also add if no overlap within 25%

def get_length(closeness_subgraph, uniprot):
    length_dict = {}
    top_percent_len = {}
    for key in closeness_subgraph:
        curr_collection= uniprot_db.query('{}'.format(key))
        #soem missing values within db so watchout for those
        for document in curr_collection:
            curr_acc = document['_id']
            curr_data = document['data']
            curr_len = curr_data['LEN']
            length_dict[curr_acc] = curr_len
            
            N = math.ceil((len(length_dict) * 0.25))
            
            top_percent_len = dict(sorted(length_dict.items(), key = lambda x: x[1], reverse = True)[:N])
            #print(top_percent_len)
            
    res = 0
    for val in length_dict.values():
        res += val
        avg = res/len(length_dict)
          
    for key, value  in top_percent_len.items():
        if value < avg:
            del top_percent_len[key]
        else:
            length_devided_byAVG = value/avg  
            #print('length', key, length_devided_byAVG)
            top_percent_len[key] = length_devided_byAVG
    return(top_percent_len )   

def potential_representative(length_dict, top_percent_cent):     
    rep_potential = []
    rep_pres = {}
    found_rep = 0
    for x in length_dict:
        for y in top_percent_cent:
            if x == y:
                rep_potential.append(x)
        
        if x in rep_potential:
            rep_pres[x] = 1
        else:
            rep_pres[x] = 0
    return(rep_pres)
    
                
                    #print("no match in centrality and length top 25%")

    #return(rep_potential)

def multiply(top_percent_cent, length_dict):
    constant_centrality = 1
    constant_length = 1
    combined_values_rep = {}
    for i in top_percent_cent:
        for j in length_dict:
            if i ==j:
                multiplied = (constant_centrality*top_percent_cent[i]) * (length_dict[j]*constant_length)
                combined_values_rep[i] = multiplied   
    return combined_values_rep

get_centrality(subgraphs)
