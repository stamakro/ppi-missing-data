import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.special import binom
from utilities import *
from sklearn.model_selection import KFold
from filereader import *

species = ['yeast', 'thalia', 'tomato']
speciesNames = ['S. cerevisiae', 'A. thaliana', 'S. lycopersicum']

metricNames = ['AUC', 'Fmax', 'Smin']
metrics = ['AUC', 'Fmax', 'Smin']

xx = np.linspace(0, 1., 51)


'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!

for i, (m, mn) in enumerate(zip(metrics, metricNames)):
	fig = plt.figure()
	for j, (s, sn) in enumerate(zip(species, speciesNames)):

		ax = fig.add_subplot(3, 1, j + 1)

		for fold in range(5):
			with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
				data = pickle.load(f)

			if fold == 0:
				expAuc = data[i][0]
	        else:
	            expAuc = np.hstack((expAuc, data[i][0]))

		ax.hist(expAuc, bins = np.linspace(0, 1, 21), density = True, color='C0', edgecolor='k')
		dens = gaussian_kde(expAuc)

		ax.plot(xx, dens(xx), color='C1')

		ax.set_xlabel(mn + ' ' + sn, fontsize=18)
		ax.set_ylabel('Probability Density', fontsize=14)

		plt.setp(ax.get_xticklabels(), fontsize=14)
		plt.setp(ax.get_yticklabels(), fontsize=14)
'''

'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!
ylims = [[0.45, 0.35, 0.40], [0.37, 0.30, 0.35], [0.68, 0.82, 1.0]]
for i, (m, mn) in enumerate(zip(metrics, metricNames)):
	fig = plt.figure()
	for j, (s, sn) in enumerate(zip(species, speciesNames)):

		if s == 'tomato':
			Nsources = 8
		else:
			Nsources = 9

		ax = fig.add_subplot(1, 3, j + 1)

		nn = np.zeros((5,))
		for fold in range(5):

			with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
				data = pickle.load(f)

			if fold == 0:
				perf = data[i]
			else:
				perf = np.hstack((perf, data[i]))

			with open('../res/' + s + '/naive/' + str(fold) + 'naive.pkl') as f:
				data = pickle.load(f)

			nn[fold] = np.mean(data[i])


		counter = 0
		meanP = np.zeros((Nsources + 1,))

		for ii in xrange(Nsources + 1):
			nreps = int(binom(Nsources, ii))

			dataC = perf[counter:counter+nreps]

			meanP[ii] = np.mean(dataC)

			ax.scatter(2*ii + 1.5 * np.random.rand(nreps) - 0.75,np.mean(dataC,1), color='C0', edgecolor='k' )
			ax.axvline(2*ii+1, color='k', linestyle=':', alpha=0.3)

			counter += nreps

		ax.plot(np.arange(0,2*(Nsources+1),2), meanP, color='C1')
		ax.set_xticks([2*k for k in range(Nsources+1)])
		ax.set_xticklabels([str(k) for k in range(Nsources+1)])

		ax.set_xlabel('Number of STRING data sources', fontsize=14)
		ax.set_ylabel(mn + ' ' + sn, fontsize=18)

		plt.setp(ax.get_xticklabels(), fontsize=11)
		plt.setp(ax.get_yticklabels(), fontsize=11)

		ax.set_ylim(0, ylims[i][j])

		ax.axhline(np.mean(nn), color='C2', linestyle='--', linewidth=1.5)
'''
'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!

xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0

xx = np.linspace(0.0, 1.0, 30)

gmd = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for i, (m, mn) in enumerate(zip(metrics, metricNames)):
	fig = plt.figure()
	for j, (s, sn) in enumerate(zip(species, speciesNames)):

		if s == 'tomato':
			ind = np.array([0, 7])
		else:
			ind = np.array([0, 5])


		for fold in range(5):

			with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
				data = pickle.load(f)

			if fold == 0:
				perf = data[i][ind]
			else:
				perf = np.hstack((perf, data[i][ind]))

		ax = fig.add_subplot(3, 2, 2 * j + 1)
		ax.scatter(perf[0], perf[1], color='C0', alpha=0.4)
		ax.plot(xx, xx, 'k--')
		ax.set_xlim(0., 1.)
		ax.set_ylim(0., 1.)
		ax.set_xlabel(gmd[2*j] + ' ' + sn + ' EXP')
		ax.set_ylabel('EXP + Text mining')


		ax = fig.add_subplot(3, 2, 2 * j + 2)
		#dd = np.array(perf)
		X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
		positions = np.vstack([X.ravel(), Y.ravel()])
		kernel = gaussian_kde(perf)
		Z = np.reshape(kernel(positions).T, X.shape)
		myplot = ax.pcolormesh(X, Y, Z, cmap='Oranges')
		fig.colorbar(myplot, ax=ax)
		ax.plot(xx, xx, 'k--')

		ax.set_xlim(0., 1.)
		ax.set_ylim(0., 1.)
		ax.set_xlabel(gmd[2*j+1] + ' ' + sn + ' EXP')
		ax.set_ylabel('EXP + Text mining')

		plt.suptitle(sn, y=-0.01)



	break

'''
'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!

fractions = np.array(range(0, 100, 10) + [95, 99, 100])

exponly = np.zeros(fractions.shape)
expString = np.zeros(fractions.shape)

ylims = [0.45, 0.45, 0.80]

for i, (m, mn) in enumerate(zip(metrics, metricNames)):
	fig = plt.figure()
	for j, frac in enumerate(fractions):
		Nsources = 9

		if frac > 0.:
			ax = fig.add_subplot(3, 4, j)

		for fold in range(5):

			with open('../res/yeast/missing' + str(int(frac)) + '/perprotein/' + str(fold) + '.pkl') as f:
				data = pickle.load(f)

			if fold == 0:
				perf = data[i]
			else:
				perf = np.hstack((perf, data[i]))


		exponly[j] = np.mean(perf[0])

		if mn == 'Smin':
			expString[j] = np.min(np.mean(perf[1:], 1))
		else:
			expString[j] = np.max(np.mean(perf[1:], 1))

		counter = 0
		meanP = np.zeros((Nsources + 1,))

		for ii in xrange(Nsources + 1):
			nreps = int(binom(Nsources, ii))

			dataC = perf[counter:counter+nreps]

			meanP[ii] = np.mean(dataC)

			if frac > 0:
				ax.scatter(2*ii + 1.5 * np.random.rand(nreps) - 0.75,np.mean(dataC,1), color='C0', edgecolor='k' )
				ax.axvline(2*ii+1, color='k', linestyle=':', alpha=0.3)

			counter += nreps

		if frac > 0.:
			ax.plot(np.arange(0,2*(Nsources+1),2), meanP, color='C1')
			ax.set_xticks([2*k for k in range(Nsources+1)])
			ax.set_xticklabels([str(k) for k in range(Nsources+1)])

			ax.set_xlabel('Number of STRING data sources', fontsize=14)
			ax.set_ylabel(mn + ' ' + str(frac), fontsize=18)

			plt.setp(ax.get_xticklabels(), fontsize=11)
			plt.setp(ax.get_yticklabels(), fontsize=11)

		#ax.set_ylim(0, ylims[i][j])

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(1,1,1)


	ax2.plot(fractions, exponly, color='C4', label='EXP')
	ax2.plot(fractions, expString, color='C9', label='EXP + FG-STRING')
	ax2.scatter(fractions, exponly, color='C4')
	ax2.scatter(fractions, expString, color='C9')
	ax2.set_ylim(0, ylims[i])

	ax2.axhline(0.37, color='k', linestyle='--', linewidth=1, label='SEQ')

	ax2.set_xlabel('% missing edges', fontsize=16)
	ax2.set_ylabel(mn, fontsize=18)

	plt.setp(ax2.get_xticklabels(), fontsize=14)
	plt.setp(ax2.get_yticklabels(), fontsize=14)

	plt.legend()
	break
sys.exit(0)
'''

'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!
errorDict = dict()
errorDict['elinewidth'] = 5
#errorDict['capsize'] = 20
errorDict['capthick'] = 2
errorDict['ecolor'] = 'k'

x = np.array([0,2,4,5,7,8])
y = np.array([0.02, 0.37, 0.42, 0.42, 0.33, 0.36])
yerr = np.array([0.001, 0.003, 0.002, 0.002, 0.002, 0.003])
cc = ['C4', 'C4', 'C0', 'C1', 'C0', 'C1']
names = ['random', 'naive', 'exp 100%', 'exp 100%\n+ STRING', 'exp 10%', 'exp 10%\n+ STRING']
names2 = [n.upper() for n in names]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.bar(x, y, width=1, yerr=yerr, color=cc, error_kw=errorDict)

ax.axhline(y[1], color='k', linestyle='--', linewidth=2)

ax.set_xticks(x)
ax.set_xticklabels(names2)

ax.set_xlabel('AFP method', fontsize=16)
ax.set_ylabel('AUC', fontsize=16)
plt.setp(ax.get_xticklabels(), fontsize=12)
plt.setp(ax.get_yticklabels(), fontsize=12)
'''
'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!

for j, (s, sn) in enumerate(zip(species, speciesNames)):
	#if j != 1:
	#	continue

	[Y, geneNames, termNames, gene2row, term2col] = goLoader(s)
	[Aexp, ppiGene2row, ppiRow2gene] = getPPInetwork(s, 'biogrid')
	assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)

	[Aexp, ppiGene2row, ppiRow2gene] = matchNetworkAndLabels(Aexp, ppiGene2row, ppiRow2gene, geneNames)

	cv = KFold(n_splits=5, shuffle=True, random_state=656391)
	np.random.seed(1901273)
	for fold, (train, test) in enumerate(cv.split(Y)):

		if fold == 0:
			degree = np.sum(Aexp[test][:, train], 1)
			assert degree.shape[0] == test.shape[0]
		else:
			degree = np.hstack((degree, np.sum(Aexp[test][:, train], 1)))


	for i, (m, mn) in enumerate(zip(metrics, metricNames)):

		fig = plt.figure()
		ax = fig.add_subplot(111)

		for fold in range(5):
			with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
				data = pickle.load(f)

			if fold == 0:
				expAuc = data[i][0]
			else:
				expAuc = np.hstack((expAuc, data[i][0]))

		X, Y = np.mgrid[0.:np.max(np.log10(degree+1)):100j, 0.:1.:100j]
		positions = np.vstack([X.ravel(), Y.ravel()])
		perf = np.vstack((np.log10(degree+1), expAuc))
		kernel = gaussian_kde(perf)
		Z = np.reshape(kernel(positions).T, X.shape)
		myplot = ax.pcolormesh(X, Y, Z, cmap='Oranges')
		fig.colorbar(myplot, ax=ax)

		ax.set_xlabel('log(degree+1) ' + sn, fontsize=18)
		ax.set_ylabel('AUC ' + sn, fontsize=14)

		plt.setp(ax.get_xticklabels(), fontsize=14)
		plt.setp(ax.get_yticklabels(), fontsize=14)
		break



'''

'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!
fig1 = plt.figure()
fig2 = plt.figure()

for si, s in enumerate(species):
	dag, mapping = read_ontology_from_file('../data/go/go-final.obo')
	[Y, geneNames, termNames, gene2row, term2col] = goLoader(s)
	[Aexp, ppiGene2row, ppiRow2gene] = getPPInetwork(s, 'biogrid')
	assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)
	[Aexp, ppiGene2row, ppiRow2gene] = matchNetworkAndLabels(Aexp, ppiGene2row, ppiRow2gene, geneNames)
	[Aexp2, ppiGene2rowS, ppiRow2geneS] = getPPInetwork(s, 'experiments')

	Aexp2 = matchNetworkAndLabels(Aexp2, ppiGene2rowS, ppiRow2geneS, geneNames)[0]
	if np.max(Aexp2) > 0.:
		threshold = np.median(Aexp2[Aexp2 > 0])
		Aexp2 = (Aexp2 >= threshold).astype(int)
		Aexp = np.maximum(Aexp, Aexp2)


	ind = np.zeros((Aexp.shape[0],), int)
	cv = KFold(n_splits=5, shuffle=True, random_state=656391)
	np.random.seed(1901273)

	c = 0
	for i, (tr, ts) in enumerate(cv.split(Y)):
		ind[c:c+ts.shape[0]] = ts
		c+= ts.shape[0]

	AA = Aexp[ind][:, ind]
	degree = np.sum(AA, 0)


	for fold in range(5):
		with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
			data = pickle.load(f)

		if fold == 0:
			expAuc = data[0][0]
		else:
			expAuc = np.hstack((expAuc, data[0][0]))


	print(s, spearmanr(degree, expAuc))


	ax = fig1.add_subplot(1,3, si+1)
	ax.scatter(np.log10(degree+1), expAuc)

	ax.set_title(s)
	ax.set_xlabel('log(degree + 1)', fontsize=16)
	ax.set_ylabel('Test AUC', fontsize=16)
	plt.setp(ax.get_xticklabels(), fontsize=12)
	plt.setp(ax.get_yticklabels(), fontsize=12)

	bins = np.array([0, 1, 5, 10, 20,1e8])

	aucs = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	for d, v in zip(degree, expAuc):
		index = 0

		while bins[index] < d:
			index += 1

		aucs[index].append(v)

	mm = [np.median(aucs[ii]) for ii in range(bins.shape[0])]
	ss = [np.std(aucs[ii], ddof=1) for ii in range(bins.shape[0])]


	ax = fig2.add_subplot(1, 3, si+1)
	ax.boxplot([aucs[k] for k in aucs])
	ax.set_xlabel('Node degree in EXP network', fontsize=16)
	ax.set_ylabel('Test AUC', fontsize=16)
	ax.set_xticklabels(['0', '1', '2-5', '6-10', '11-20', '>20'])
	ax.set_title(s)

	plt.setp(ax.get_xticklabels(), fontsize=12)
	plt.setp(ax.get_yticklabels(), fontsize=12)


'''
#'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!
fig1 = plt.figure()
fig2 = plt.figure()

indexGmd = [4, 4, 3]

aa = ['(a) ', '(b) ', '(c) '] 

for si, s in enumerate(species):
	dag, mapping = read_ontology_from_file('../data/go/go-final.obo')
	[Y, geneNames, termNames, gene2row, term2col] = goLoader(s)
	[Aexp, ppiGene2row, ppiRow2gene] = getPPInetwork(s, 'biogrid')
	assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)
	[Aexp, ppiGene2row, ppiRow2gene] = matchNetworkAndLabels(Aexp, ppiGene2row, ppiRow2gene, geneNames)
	[Aexp2, ppiGene2rowS, ppiRow2geneS] = getPPInetwork(s, 'experiments')

	Aexp2 = matchNetworkAndLabels(Aexp2, ppiGene2rowS, ppiRow2geneS, geneNames)[0]
	if np.max(Aexp2) > 0.:
		threshold = np.median(Aexp2[Aexp2 > 0])
		Aexp2 = (Aexp2 >= threshold).astype(int)
		Aexp = np.maximum(Aexp, Aexp2)


	ind = np.zeros((Aexp.shape[0],), int)
	cv = KFold(n_splits=5, shuffle=True, random_state=656391)
	np.random.seed(1901273)

	c = 0
	for i, (tr, ts) in enumerate(cv.split(Y)):
		ind[c:c+ts.shape[0]] = ts
		c+= ts.shape[0]

	AA = Aexp[ind][:, ind]
	degree = np.sum(AA, 0)


	for fold in range(5):
		with open('../res/' + s + '/missing0/perprotein/' + str(fold) + '.pkl') as f:
			data = pickle.load(f)

		if fold == 0:
			expAuc = data[0][indexGmd[si]] - data[0][0]
		else:
			expAuc = np.hstack((expAuc, data[0][indexGmd[si]]-data[0][0]))


	print(s, spearmanr(degree, expAuc))


	ax = fig1.add_subplot(1,3, si+1)
	ax.scatter(np.log10(degree+1), expAuc)

	ax.set_title(s)
	ax.set_xlabel('log(degree + 1)', fontsize=16)
	ax.set_ylabel('Test AUC', fontsize=16)
	plt.setp(ax.get_xticklabels(), fontsize=12)
	plt.setp(ax.get_yticklabels(), fontsize=12)

	bins = np.array([0, 1, 5, 10, 20,1e8])

	aucs = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	for d, v in zip(degree, expAuc):
		index = 0

		while bins[index] < d:
			index += 1

		aucs[index].append(v)

	mm = [np.median(aucs[ii]) for ii in range(bins.shape[0])]
	ss = [np.std(aucs[ii], ddof=1) for ii in range(bins.shape[0])]


	ax = fig2.add_subplot(1, 3, si+1)
	ax.boxplot([aucs[k] for k in aucs])
	ax.set_xlabel(aa[si] + 'Node degree in EXP network', fontsize=16)
	ax.set_ylabel('Test AUC', fontsize=16)
	ax.set_xticklabels(['0', '1', '2-5', '6-10', '11-20', '>20'])
	ax.set_title(speciesNames[si],fontsize=18)

	plt.setp(ax.get_xticklabels(), fontsize=14)
	plt.setp(ax.get_yticklabels(), fontsize=14)

#'''
'''
#!!!!!!!!!!!!!!!!!!!!!  FIGURE  !!!!!!!!!!!!!!!!!!!!!!!
data = np.array([[0.02, 0.37, 0.42, 0.42, 0.37, 0.38], [0.01, 0.33, 0.16, 0.30, 0.33, 0.34], [0.10, 0.39, 0.11, 0.3, 0.42, 0.42]])
err = np.array([[0.001, 0.003, 0.002, 0.002, 00., 0.003], [0.001, .002, 0.003, 0.001, 0.002, 0.002], [0.005, 0.009, 0.004, 0.006, 0.008, 0.008]])
methods = ['random', 'naive', 'EXP', 'EXP+FG-STRING','SEQ', 'EXP+SEQ']


fig = plt.figure()
ax = fig.add_subplot(111)
colors = ['k', 'k', 'lightgray', 'dimgray', 'w', 'w']
hatches = ['//', None, '\\\\', '\\\\', None, '\\\\']

ec = ['w', 'k', 'k', 'k', 'k', 'k']


for i in range(6):
	ax.bar([7*j+i for j in range(3)], data[:,i], yerr=err[:,i], color=colors[i], hatch=hatches[i], label= methods[i], edgecolor=ec[i])
	ax.set_xticks([2.5+7*j for j in range(3)])
	ax.set_xticklabels(speciesNames)

	ax.set_ylabel('AUC', fontsize=16)
	ax.set_ylabel('AUC', fontsize=16)

	plt.setp(ax.get_xticklabels(), fontsize=16)
	plt.setp(ax.get_yticklabels(), fontsize=13)

plt.legend(loc='upper center')
'''






plt.show()
