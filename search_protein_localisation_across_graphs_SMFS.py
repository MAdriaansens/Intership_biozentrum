import os
import glob
import sys
import subprocess

infolder = sys.argv[1]
in_protein = sys.argv[2]


directionairy_list = list(glob.glob('{}/*'.format(infolder)))
community_dir_id_list = []

for i in directionairy_list:
    community_inof = '{}/community_id.txt'.format(i)
    community_dir_id_list.append(community_inof)

word = "{}".format(in_protein)


for community_txt in community_dir_id_list:
    cluster =community_txt.split("/")[-2].split("/")[-1]
    with open('{}'.format(community_txt), 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
        for line in lines:
            # check if string present on a current lin
            if line.find(word) != -1:
                print(word, 'string exists in file')
                #print('Line Number:', lines.index(line))
                line =  line.split('{')[0]
                print('Community_info for {}'.format(cluster), line)
