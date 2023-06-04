import sys
sys.path.append('/scicore/home/schwede/soares0000/projects/dark_protein_universe/my_menzi/')
 
from src import extract_uniprot  as uniprot
 
MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077
 
uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
