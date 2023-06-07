#make graph to see localisation of different connected components
import datashader as ds, pandas as pd, colorcet as cc

infile = pd.read_csv('Archaea_comm_edge_list_si50_cov80-metadata.csv', sep = ';')
list_cc = []

infile['id'] = infile['id'].astype(str) 
infile['x'] = infile['x'].astype(float)
infile['y'] = infile['y'].astype(float)

for i in infile['id']:
    i = i.split('[')[0].split(']')[-1]
    list_cc.append(i)
    


df = pd.DataFrame(list_cc)

infile['cc'] = df
    #print(infile['cc'])
    #infile = infile[(infile['cc']== '0') |(infile['cc']== '1') |(infile['cc']== '2') |(infile['cc']== '3')|(infile['cc']== '4')]
    #print(infile)
    
agg = ds.Canvas().points(infile, 'x', 'y')
ds.tf.set_background(ds.tf.shade(agg, cmap=cc.isolum), "black")
