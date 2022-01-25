from Bio import Entrez
import time
import pandas as pd

Entrez.email = ""

def search_literature(term):
	store_papers = []
	handle = Entrez.esearch(db='pubmed', term=term,retmax=10000)
	pubmed_out = Entrez.read(handle)
	pmid = ",".join(map(str,pubmed_out['IdList']))
	handle = Entrez.esummary(db='pubmed', id=pmid,retmax=10000)
	paper_deets = Entrez.read(handle)
	for paper in paper_deets:
		store_papers.append(paper['Title'] + "\t" + paper['LastAuthor'] + "\t" + paper['PubDate'] + "\t" + paper['FullJournalName'] + "\t" + paper['Id'])
	return store_papers


def write_to_file(filename, listy):
	with open(filename, mode="w") as txtfile:
		for item in listy:
			txtfile.write("%s \n" % item)

def main():
	gene_term = 'ASCC1[Title] | ASCC1 muscle'
	disorder_term = 'centronuclear[Title] myopathy[Title]'
	gene_phenotype_term = 'ASCC1 myopathy'

	gene_papers = search_literature(gene_term)
	disorder_papers = search_literature(disorder_term)
	gene_phen_papers = search_literature(gene_phenotype_term)
	write_to_file("disease_papers.txt", disorder_papers)
	write_to_file("ascc1_papers.txt", gene_papers)
	write_to_file("ascc1_myopathy_papers.txt", gene_phen_papers)

if __name__ == "__main__":
	main()
