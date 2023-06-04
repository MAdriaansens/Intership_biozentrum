import os
import pandas as pd
import sys
sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')

MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077
from src import extract_uniprot  as uniprot
uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)


uniprot_list = 'databases/uniprot-compressed_true_download_true_fields_id_2Cidentity_format_ts-2023.06.01-15.11.31.71.tsv'

def get_sequences_uniprot(uniprot_list):
    
    #get uniprot id into df
    df = pd.read_table(uniprot_list)
    uniprot_id = []
    for i in df['Cluster ID']:
        uniprot_id.append(i.split('_')[-1])
    df_uniprot_id = pd.DataFrame(uniprot_id, columns = ['uniprot_id'])
    #print(df_uniprot_id)
    
    df_uniprot_id.to_csv('uniprot_ids.csv', index = False)

get_sequences_uniprot(uniprot_list)

df_uniprot_id = pd.read_table('uniprot_ids.csv')

def search_uniprot(uniprot_db, df_uniprot_id):
    #print(df_uniprot_id)
    #search against uniprot for sequences
    print("started")
    count = 0
    for i in df_uniprot_id['uniprot_id']:
        curr_collection= uniprot_db.query('{}'.format(i))
        for document in curr_collection:
            #print(document
            curr_acc = document['_id']
            curr_data = document['data']
            curr_taxid = curr_data['TAXID']
            curr_name = curr_data['NAME']
            curr_title = curr_name['TITLE']
            curr_seq = curr_data['SEQ']
            curr_len = curr_data['LEN']
            if curr_len > 14:
                outfile = open('Archaea_uniref50.fasta', 'a')
                outfile.write(">"+"{}".format(curr_acc) + "{}".format(curr_name) + '\n' + '{}'.format(curr_seq) + "\n")
                outfile.close()
                count = count +1
                if count % 10000 == 0:
                    percentage = (count/1939015 * 100)
                    print("we are at protein {}/1939015 which is {} %".format(count, percentage))
    print("finished")

search_uniprot(uniprot_db, df_uniprot_id)
