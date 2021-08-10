from Bio import SeqIO
import datetime as dt
import csv
from collections import defaultdict

def decimal_date(date_string):

    date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    final_date = date.year + float(date.toordinal() - start) / year_length

    return str(final_date)

# def get_options_per_tree(trait_locs, trait_dict, taxon_set_file):

#     options_per_tree = defaultdict(dict)

#     seq_to_set = {}
#     with open(taxon_set_file) as f:
#         data = csv.DictReader(f)
#         for line in data:
#             seq = line['sequence_name']
#             taxon_set = line["taxon_set"]

#     for trait, position in trait_locs.items():


def process_coordinates(trait_file):

    trait_dict = defaultdict(dict)

    with open(trait_file) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            inner_dict = {}
            inner_dict["longitude"] = l['longitude']
            inner_dict["latitude"] = l["latitude"]
            inner_dict["coordinates"] = f"{l['latitude']} {l['longitude']}"

            trait_dict[l['taxon']] = inner_dict

    return trait_dict
            
             








        

    




