#make graph to see localisation of different connected components
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

infile = pd.read_csv('results/Archaea_comm_edge_list_si50_cov80-metadata.csv', sep = ';')
list_cc = []

#entire panda columns are string or float however to assure correct editing astype() is used 
infile['id'] = infile['id'].astype(str) 
infile['x'] = infile['x'].astype(float)
infile['y'] = infile['y'].astype(float)

for i in infile['id']:
    i = i.split('[')[0].split(']')[-1]
    list_cc.append(i)
    


df_cc = pd.DataFrame(list_cc)

infile['cc'] = df_cc
#specify of which connected component you wish a plot, is the line below is excluded a plot of the entire graph will be made
infile = infile[(infile['cc']== '50')]
sns.scatterplot(x = infile['x'], y = infile['y'], hue= infile['cc'])
plt.show()
#print('if this shows without plot you have been scammed')
