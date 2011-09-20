#cobra.io.resolver.py
#Parsing files from the Rosetta Resolver.
from collections import defaultdict
print "WARNING: Functions in cobra.io.resolver are not yet functional."
def parse_entrez(file_name):
    """Parses a file produced from Rosetta resolver where the expression values are
    collapsed to the entrez_gene_ids.

    TODO: Create the nomenclature tables from the slide annotation.
    """
    field_to_index = {'log_ratio': 5,
                      'log_error': 6,
                      'p_value': 7,
                      'intensity_1': 8,
                      'intensity_2': 9}
    with open(file_name) as in_file:
        in_file.readline() #Skip the header row.
        the_data = [x.strip().split('\t')
                    for x in in_file.readlines()]
        while len(the_data[-1]) < 10:
            the_data.pop()
    entrez_gene_id_index = 3
    data_dict = defaultdict(dict)
    #assume that entrez_gene_ids are unique in the files.
    for the_row in the_data:
        tmp_dict = data_dict[int(the_row[entrez_gene_id_index])]
        [tmp_dict.update({k: float(the_row[v])})
         for k, v in field_to_index.items()]
    return(data_dict)

