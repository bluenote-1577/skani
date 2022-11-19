import numpy as np
from collections import defaultdict
import seaborn
import matplotlib.pyplot as plt
import sys
sys.setrecursionlimit(100000)
from scipy.cluster import hierarchy
import scipy
file = sys.argv[1]
counter = 0
items = 0
labels = []
condensed = []
matrix = []
all_labels = set()
delim = '\t'
#delim = ','
for line in open(file, 'r'):
    if counter == 0:
        #print(line)
        spl = line.split(delim)
        if len(spl) > 2:
            items = len(spl)
        else:
            items = int(line.split(delim)[-1])
        print(items)
        matrix = [[] for x in range(items)]
        counter += 1
        continue
    if delim in line:
        spl = line.split(delim);
    else:
        spl = line.split();
    print(spl[0].split('/')[-1])
    labels.append(spl[0].split('/')[-1])
    endpoints = range(1,counter)
    for i in endpoints:
        matrix[i-1].append(float(spl[i]))
    counter += 1

for vec in matrix:
    for score in vec:
        condensed.append(1 - score)


cmap = seaborn.cm.rocket_r

#Z = hierarchy.linkage(condensed, 'single')
#Z = hierarchy.linkage(condensed, 'complete')
Z = hierarchy.linkage(condensed, 'average')
square_mat = scipy.spatial.distance.squareform(condensed)
if len(sys.argv) > 2:
    vmax = float(sys.argv[2])
    cg = seaborn.clustermap(square_mat, row_linkage = Z, col_linkage = Z, vmax = vmax, cmap = cmap)
else:
    cg = seaborn.clustermap(square_mat, row_linkage = Z, col_linkage = Z, cmap = cmap)

 
print(cg.dendrogram_row.reordered_ind)
re = [labels[x] for x in cg.dendrogram_row.reordered_ind]
xticks = [x for x in range(len(labels))]
cg.ax_heatmap.set_xticks(xticks)
cg.ax_heatmap.set_xticklabels(re, rotation=90)
#cg.ax_heatmap.set_yticklabels(re, rotation=0)
plt.tight_layout()
plt.show()
