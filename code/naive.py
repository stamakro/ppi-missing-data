import numpy as np
import pickle
import sys
from utilities import *
from filereader import read_ontology_from_file
from gyros import getParentsCoord, calculateIC

species = sys.argv[1]


dag, mapping = read_ontology_from_file('../data/go/go-final.obo')

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)


cv = KFold(n_splits=5, shuffle=True, random_state=656391)

np.random.seed(1901273)

thresholds = np.linspace(0, 1.0, 21)

for fold, (train, test) in enumerate(cv.split(Y)):

	print fold

	termNames = np.array(termNames)

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

	freqs = np.sum(Ytrain, 0) / float(Ytrain.shape[0])
	assert freqs.shape[0] == Ytest.shape[1]

	Y_posteriors = np.tile(freqs, (Ytest.shape[0],1))

	pc_auprc = average_precision_score(Ytest.T, Y_posteriors.T, average=None)

	pc_f1 = np.zeros((thresholds.shape[0], Ytest.shape[0]))
	pc_nsd = np.zeros((thresholds.shape[0], Ytest.shape[0]))

	for i, thres in enumerate(thresholds):
		Ypred = (Y_posteriors >= thres).astype(int)

		pc_nsd[i] = normalizedSemanticDistance(Ytest, Ypred, ic)[2]
		pc_f1[i] = f1_score(Ytest.T, Ypred.T, average=None)


	with open('../res/' + species + '/naive/' + str(fold) + 'naive.pkl', 'w') as f:
		pickle.dump([pc_auprc, pc_f1[np.argmax(np.mean(pc_f1,1))], pc_nsd[np.argmin(np.mean(pc_nsd,1))]], f)
