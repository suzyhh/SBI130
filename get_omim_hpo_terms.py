from Bio import Entrez
import phizz, time
import pandas as pd
from hpo import HPOApi
import seaborn as sns
import matplotlib.pyplot as plt

Entrez.email = ""

def get_omim_ids(database, search_term):
	'''
	Function to obtain omim disease IDs associated with given terms
	'''
	handle = Entrez.esearch(db=database, term=search_term, retmax=10000)
	record = Entrez.read(handle)

	omim_str = ",".join(map(str,record['IdList']))
	
	handle = Entrez.esummary(db="omim", id=omim_str)
	omim_details  = Entrez.read(handle)

	disease_code = {}

	for i in omim_details:
		disease_code[i["Id"]] = i["Title"]
		# save this as a tsv file

	omim_df = pd.DataFrame.from_dict(disease_code, orient="index").reset_index()
	omim_df.columns = ['omim', 'description']
	omim_df.to_csv('omim_ids.csv', sep ='\t', index=False)

	omim_ids = record["IdList"]
	return omim_ids

def get_hpo_terms(omim_ids, min_shared):
	'''
	Get the HPO terms associated with the OMIM gene IDs
	'''
	hpo_dict = {'hpo': [], 'description': []}
	hpo_terms = []
	
	for omim in omim_ids:
		ind_hpo_terms = []
		search_term = "OMIM:" + omim
		phizz_dict = phizz.query_disease([search_term])
		for hpo in phizz_dict:
			if hpo['hpo_term'] not in ind_hpo_terms:
				hpo_terms.append(hpo['hpo_term'])
				ind_hpo_terms.append(hpo['hpo_term'])

	# To filter the list, count the occurances of each HPO and only keep those that appear in >3 OMIM diseases
	hpo_counts = dict()
	for i in hpo_terms:
		hpo_counts[i] = hpo_counts.get(i, 0) + 1

	shared_terms = list({x: count for x, count in hpo_counts.items() if count >=min_shared}.keys())
	common = ['HP:0000006', 'HP:0000007']
	shared_terms = [x for x in shared_terms if x not in common]


	for term in shared_terms:
		query_dict = phizz.query_hpo([term])
		for res in query_dict:
			hpo_dict['hpo'].append(res['hpo_term'])
			hpo_dict['description'].append(res['description'])
	
	shared_hpo_df = pd.DataFrame(data=hpo_dict)
	shared_hpo_df.columns = ['hpo', 'description']
	shared_hpo_df.to_csv('disease_hpo_terms.csv', sep =',', index=False)

	return shared_terms

def get_gene_hpo_terms(suid_file):
	'''
	Get the HPO terms associated with each gene identified by Moduland
	'''
	hpo = HPOApi()

	with open(suid_file) as f:
		gene_ids = [line.strip() for line in f.readlines()[1:]]
	
	gene_terms = []

	for id in gene_ids:
		try:
			gene_terms.append(hpo.gene_search(id))
		except:
			pass
	
	return gene_terms

def create_matrixes(gene_terms, disease_terms):
	date_today = time.strftime('%Y%m%d')
	to_plot = {gene.gene_symbol: [term.ontology_id for term in gene.terms] for gene in gene_terms}

	values = list(set([ x for y in to_plot.values() for x in y]))
	data = {}
	for key in to_plot.keys():
		data[key] = [ 1 if value in to_plot[key] else 0 for value in values ]

	df = pd.DataFrame(data, index=values).transpose()
	df.to_csv(date_today + 'df.csv', sep =',')
	incidence_df = df.filter(items=disease_terms)
	incidence_df = incidence_df.loc[~(incidence_df==0).all(axis=1)]
	incidence_df.to_csv(date_today + '_incidence.csv', sep =',')

	plt.rcParams.update({'font.size': 18})
	fig, ax = plt.subplots(figsize=(30,40))
	sns.heatmap(incidence_df, cmap="viridis")
	plt.savefig(date_today + "_incidence.png")

	adjacency_df = incidence_df.dot(incidence_df.transpose())
	fig, ax = plt.subplots(figsize=(50,50))
	sns.heatmap(adjacency_df, cmap="viridis")
	plt.savefig(date_today + "_adjacency.png")
	adjacency_df.to_csv(date_today + '_adjacency_matrix.csv', sep =',')

	plt.rcParams.update({'font.size': 8})
	fig, ax = plt.subplots(figsize=(50,50))
	sns.clustermap(adjacency_df, cmap="viridis", xticklabels=True, yticklabels=True)
	plt.savefig(date_today + '_adjacency_clustermap.png')

	return adjacency_df

def filter_adjacency_matrix(adjacency_df):
	date_today = time.strftime('%Y%m%d')
	panel_genes = ['ACTA1','ACTN2','ADSSL1','BIN1','CACNA1S','CCDC78','CFL2','COL12A1','COL6A1','COL6A2','COL6A3','DNM2','DOK7','ECEL1','EPG5','FKBP14','FXR1','KBTBD13','KLHL40','KLHL41','LMNA','LMOD3','MAP3K20','MEGF10','MICU1','MTM1','MYBPC1','MYH2','MYH3','MYH7','MYH8','MYL1','MYMK','MYO18B','MYPN','NEB','ORAI1','PAX7','PIEZO2','PYROXD1','RYR1','RYR3','SCN4A','SELENON','SLC25A4','SPEG','STAC3','STIM1','TNNI2','TNNT1','TNNT3','TPM2','TPM3','TRIP4','TTN','VMA21']
	temp_filtered = adjacency_df[adjacency_df.columns.intersection(panel_genes)]
	filtered_adj = temp_filtered[~temp_filtered.index.isin(panel_genes)]
	filtered_adj.to_csv(date_today + '_filtered_adjacency_matrix.csv', sep =',')
	
	plt.rcParams.update({'font.size': 10})
	fig, ax = plt.subplots(figsize=(50,50))
	sns.clustermap(filtered_adj, annot=True, cmap="viridis", xticklabels=True, yticklabels=True)
	plt.savefig(date_today + '_filtered_adjacency_clustermap.png')

	novel_sum = filtered_adj.sum(axis = 1)
	novel_sum = novel_sum.sort_values(ascending=True)
	novel_sum.to_csv(date_today + '_terms_per_gene.csv', sep=',')

	plt.rcParams.update({'font.size': 22})
	fig, ax = plt.subplots(figsize=(50,20))
	novel_sum.plot.bar()
	plt.savefig(date_today + '_terms_per_gene.png')

def main():
	omim_ids = get_omim_ids('omim', 'congenital[Title] | nemaline[Title] myopathy[Title]')
	disease_terms = get_hpo_terms(omim_ids, 4)
	gene_terms = get_gene_hpo_terms("suid_output.csv")
	adjacency_df = create_matrixes(gene_terms, disease_terms)
	filter_adjacency_matrix(adjacency_df)


if __name__ == "__main__":
	main()
