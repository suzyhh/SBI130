# SBI130

Scripts written to fulfill the SBI130 competencies.

1. get_omim_hpo_terms.py
  * Obtains a list of OMIM IDs linked to a specified search term
  * Uses the OMIM IDs to obtain a list of linked disease HPO terms
  * Takes a CSV file of Ensembl gene IDs (obtained from Cytoscape network analysis) and obtains a list of HPO terms linked to each
  * Intersects the gene and disease HPO terms
  * Builds incidence and adjacency matrices to assess the strengh of relationship between each gene and congenital myopathy and also between PanelApp genes and novel candidate genes

2. search_snp_db.py
  * Searches ClinVar for variants in a specified gene (ASCC1) and extracts relevant information
  * dbSNP is searched to obtain gnomAD allele frequencies for variants with RS IDs. 
  * Writes details to a text file

3. search_literature.py
  * Uses a specified searchterm to query PubMed for matching literature. Writes details to a text file
