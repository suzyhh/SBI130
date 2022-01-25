from Bio import Entrez
import csv

Entrez.email = ""


def search_clinvar(search_term):
	'''
	Search the clinvar database for any variants in the gene specified by search_term.
	Extract variants with relvant metadata into a dictionary
	'''

	handle = Entrez.esearch(db='clinvar', term=search_term, retmax=10000)
	record = Entrez.read(handle)
	snp_ids = ",".join(map(str,record['IdList']))

	handle = Entrez.esummary(db='clinvar', id=snp_ids)
	record2 = Entrez.read(handle)

	clinvar_dict = {}

	# iterate through the horrible Entrez dictionary object								
	for k,v in record2.items():
		for k1,v1 in v.items():
			if type(v1) == list:
				for dict in v1:
					hgvs = dict['title']
					var_type = dict['obj_type']
					path = dict['clinical_significance']['description']
					for trait_dict in dict['trait_set']:
						trait = trait_dict['trait_name']
					for dict2 in dict['variation_set']:
						if len(dict2['variation_xrefs']) > 0:
							for db_source in dict2['variation_xrefs']:
								if db_source['db_source'] == 'dbSNP':
									rs = db_source['db_id']
								else:
									rs = 'NA'
						else:
							rs = 'NA'
						clinvar_dict.update({hgvs: {
							'var_type': var_type,
							'pathogenicity': path,
							'disorder': trait,
							'rs_id': rs,
							'gnomad_freq': 'NA'
							}})
	return clinvar_dict

def search_dbsnp(dictionary):
	'''
	For the variants that have rs IDs, search the dbSNP database
	Update the dictionary to contain gnomad_exome MAF
	'''
	# get a list of rs_ids from variants with a dbSNP entry
	rsid = []
	for k,v in dictionary.items():
		if v['rs_id'] != 'NA':
			rsid.append(v['rs_id'])
	rs_str = ",".join(map(str,rsid))

	# search this list through dbsnp to obtain allele freq
	handle = Entrez.esummary(db='snp', id=rs_str)
	record3=Entrez.read(handle)

	gnomad_dict = {}

	# iterate through the slightly-less-horrible Entrez dictionary object
	for k,v in record3.items():
		for k1,v1 in v.items():
			if type(v1) == list:
				for dict in v1:
					snpid = dict['SNP_ID']
					gnomad_dict[snpid] = 'not reported' # some SNPs will have an rs ID but won't be reported in gnomad_exomes
					for i in dict['GLOBAL_MAFS']:
						if i['STUDY'] == 'GnomAD_exomes': # if the SNP has MAF data AND has a gnomad_exome MAF, update freq in dictionary
							temp_f = i['FREQ'].replace('/','=')
							freq = temp_f.split('=')[1]
							gnomad_dict[snpid] = freq

	# add allele freq to clinvar_dict
	for k,v in dictionary.items():		
		if v['rs_id'] in gnomad_dict:
			v['gnomad_freq'] = gnomad_dict[v['rs_id']]

	# add allele freq to clinvar_dict
	return dictionary

def write_to_csv(outname, dict):
	'''
	Write the dictionary to a CSV file
	'''

	# define headers
	fnames = ['hgvs','var_type','pathogenicity','disorder','rs_id', 'gnomad_freq']
	
	with open(outname, 'w') as f:
		writer = csv.DictWriter(f, fieldnames=fnames)
		writer.writeheader()
		for k in dict:
			writer.writerow({field: dict[k].get(field) or k for field in fnames})

def main():
	search_term = 'ASCC1[Gene]'
	clinvar_dict = search_clinvar(search_term)
	snp_dict = search_dbsnp(clinvar_dict)
	write_to_csv("ASCC1_SNPs.csv", snp_dict)

if __name__ == "__main__":
	main()
