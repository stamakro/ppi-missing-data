import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys

species = 'tomato'

if species == 'tomato':
    datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])
else:
    datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

ds2 = []
for d in datasources:
    ff = d.split('_')
    if len(ff) == 1:
        ds2.append(d)
    else:
        ds2.append(ff[0] + '\n' + ff[1])

datasources = ds2

Nds = len(datasources)

if species == 'yeast':
    missing = '90'
else:
    missing = '0'

fold = 0
with open('../res/' + species + '/missing' + missing + '/perprotein/' + str(fold) + '.pkl') as f:
    data = pickle.load(f)

auc = np.mean(data[0],1)

d1 = auc[1:1+Nds]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.bar(np.arange(Nds), d1)

ax.bar(-1, auc[0], color='C1')
ax.bar(Nds, np.max(auc), color='C1')
ax.set_xticks(np.arange(-1, 1+Nds))
ax.set_xticklabels(['EXP'] + datasources + ['best'])


plt.show()
