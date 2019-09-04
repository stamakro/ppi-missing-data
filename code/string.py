import numpy as np
import pickle
import sys
from utilities import *
from itertools import combinations
from sklearn.metrics import average_precision_score
from filereader import read_ontology_from_file
from gyros import getParentsCoord, calculateIC
from scipy.special import binom


species = sys.argv[1]

if len(sys.argv) == 2:
	missingFraction = 0.0
else:
	missingFraction = float(sys.argv[2])


dag, mapping = read_ontology_from_file('../data/go/go-final.obo')

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)

print('Loading experimental PPI network...')

[Aexp, ppiGene2row, ppiRow2gene] = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)

[Aexp, ppiGene2row, ppiRow2gene] = matchNetworkAndLabels(Aexp, ppiGene2row, ppiRow2gene, geneNames)

[Aexp2, ppiGene2rowS, ppiRow2geneS] = getPPInetwork(species, 'experiments')
assert np.max(np.abs((Aexp2 - Aexp2.T)) < 1e-10)
Aexp2 = matchNetworkAndLabels(Aexp2, ppiGene2rowS, ppiRow2geneS, geneNames)[0]

if np.max(Aexp2) > 0.:

	threshold = np.median(Aexp2[Aexp2 > 0])
	Aexp2 = (Aexp2 >= threshold).astype(int)
	Aexp = np.maximum(Aexp, Aexp2)



if species != 'tomato':

	datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

else:

	datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])
#assert datasource.shape[0] == 9




with open('../data/' + species + '/interactions/biogrid-final/row2protein.pkl') as f:
	ppiRow2gene = pickle.load(f)

ppiGene2row = dict()
for k in ppiRow2gene:
	ppiGene2row[ppiRow2gene[k]] = k

for i, ds in enumerate(datasources):
	print 'Reading predicted', i
	with open('../data/' + species + '/interactions/string-final/' + ds + '.pkl') as f:
		At = pickle.load(f).toarray()
		At += At.T

		assert np.max(np.abs((At - At.T)) < 1e-10)

		At = matchNetworkAndLabels(At, ppiGene2row, ppiRow2gene, geneNames)[0]

	if i == 0:
		AA = np.zeros((datasources.shape[0], At.shape[0], At.shape[1]))

	AA[i] = At

print 'Correcting for homology...'

#homology correction
homInd = np.where(datasources == 'homology')[0]

tmInd = np.where(datasources == 'textmining')[0]

AA[tmInd] = AA[tmInd] * (1 - AA[homInd])

coInd = np.where(datasources == 'cooccurence')[0]
AA[coInd] = AA[coInd] * (1 - AA[homInd])



Aexp = removeEdges(Aexp, missingFraction)


F = len(datasources)
ss = 0
for i in range(1, F + 1):
	ss += binom(F, i)

ss = int(ss)

n_folds = 5


cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

termNames = np.array(termNames)
print type(termNames)
for fold, (train, test) in enumerate(cv.split(Y)):

	Ytrain = Y[train]
	Ytest = Y[test]

	print type(termNames)
	nonempty = np.where(np.sum(Ytest, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]
	print type(termNames)
	termNames = termNames[nonempty]

	nonempty = np.where(np.sum(Ytrain, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]
	termNames = list(termNames[nonempty])

	auc = np.zeros((ss + 1, Ytest.shape[0]))
	nsmin = np.zeros(auc.shape)
	fmax = np.zeros(auc.shape)

	parentsCoord = getParentsCoord(termNames, 'P', dag, mapping)
	ic = calculateIC(Ytrain, parentsCoord)


	counter = 0
	for i in range(F + 1):

		if i == 0:

			auc[counter], fmax[counter], nsmin[counter] = predict(Aexp[test][:, train], Ytrain, Ytest, ic)
			break

		else:

			for j,p in enumerate(combinations(range(F), i)):
				counter += 1
				print counter
				if counter < start:
					continue

				if counter > start + 99:

					with open('../res/missing' + str(int(100 * missingFraction)) + '/perprotein/'+ str(start) + '.pkl', 'w') as fw:
						pickle.dump([auc, sd], fw)

					sys.exit(0)

				Apred = integrateStringScores(AA[list(p)])

				threshold = np.median(Apred[Apred > 0])
				Apred = (Apred >= threshold).astype(int)

				A = np.maximum(Aexp, Apred)


				auc[counter], fmax[counter], nsmin[counter] = predict(A[test][:, train], Ytrain, Ytest, ic)

	break

#with open('../res/missing' + str(int(100 * missingFraction))  + '/perprotein/'+ str(start) + '.pkl', 'w') as fw:
#	pickle.dump([auc, sd], fw)
