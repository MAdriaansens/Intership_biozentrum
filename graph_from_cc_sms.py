#make graph to see localisation of different connected components
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

infile = pd.read_csv('results/Archaea_comm_edge_list_si50_cov80-metadata.csv', sep = ';')
list_cc = []
infile['id'] = infile['id'].astype(str) 
infile['x'] = infile['x'].astype(float)
infile['y'] = infile['y'].astype(float)

for i in infile['id']:
    i = i.split('[')[-1].split(']')[0]
    list_cc.append(i)
    


df_cc = pd.DataFrame(list_cc)

infile['cc'] = df_cc
infile = infile[(infile['cc']== '50')]
print(infile)
sns.scatterplot(x = infile['x'], y = infile['y'], hue= infile['cc'])
plt.show()
#print('if this shows without plot you have been scammed')
