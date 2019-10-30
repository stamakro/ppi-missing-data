import numpy as np
import pickle
import sys
from utilities import *
from itertools import combinations
from sklearn.metrics import average_precision_score, roc_auc_score
from filereader import read_ontology_from_file
from gyros import getParentsCoord, calculateIC, normalizedSemanticDistance
from scipy.special import binom

species = sys.argv[1]

dag, mapping = read_ontology_from_file('../data/go/go-final.obo')

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)


print('Loading experimental PPI network...')
[Aexp, ppiGene2row, ppiRow2gene] = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)

[Aexp, ppiGene2row, ppiRow2gene] = matchNetworkAndLabels(Aexp, ppiGene2row, ppiRow2gene, geneNames)

missingFraction = 0.0
Aexp = removeEdges(Aexp, missingFraction)


[Aexp2, ppiGene2rowS, ppiRow2geneS] = getPPInetwork(species, 'experiments')
assert np.max(np.abs((Aexp2 - Aexp2.T)) < 1e-10)
Aexp2 = matchNetworkAndLabels(Aexp2, ppiGene2rowS, ppiRow2geneS, geneNames)[0]

if np.max(Aexp2) > 0.:

	threshold = np.median(Aexp2[Aexp2 > 0])
	Aexp2 = (Aexp2 >= threshold).astype(int)
	Aexp = np.maximum(Aexp, Aexp2)

with open('../data/' + species + '/interactions/biogrid-final/row2protein.pkl') as f:
	ppiRow2gene = pickle.load(f)

ppiGene2row = dict()
for k in ppiRow2gene:
	ppiGene2row[ppiRow2gene[k]] = k


Adl = getPPInetwork(species, 'dl')[0]
assert Adl.shape == Aexp.shape

Aexpdl = np.maximum(Aexp, Adl)


cv = KFold(n_splits=5, shuffle=True, random_state=656391)

np.random.seed(1901273)

termNames = np.array(termNames)
for fold, (train, test) in enumerate(cv.split(Y)):
	print fold

	Ytrain = Y[train]
	Ytest = Y[test]

	nonempty = np.where(np.sum(Ytest, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]
	termNames2 = termNames[nonempty]

	nonempty = np.where(np.sum(Ytrain, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]
	termNames2 = list(termNames2[nonempty])


	parentsCoord = getParentsCoord(termNames2, 'P', dag, mapping)
	ic = calculateIC(Ytrain, parentsCoord)

	stuff = predict(Adl[test][:, train], Ytrain, Ytest, ic)


	with open('../res/' + species + '/dl/' + str(fold) + 'dlOnly.pkl', 'w') as f:
		pickle.dump(stuff, f)

	stuff = predict(Aexpdl[test][:, train], Ytrain, Ytest, ic)

	with open('../res/' + species + '/dl/' + str(fold) + 'dlExp.pkl', 'w') as f:
		pickle.dump(stuff, f)
