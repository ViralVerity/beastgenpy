from Bio import SeqIO
from collections import defaultdict



#This is not complete - not sure how this is going to come out of the DTA
#but basically what needs to come out of this is the taxon_sets dict of name:seqs
def parse_sets(fasta_file, taxon_set_file):

    aln = SeqIO.parse(fasta_file, 'fasta')

    sequence_names = []

    for seq in aln:
        sequence_names.append(seq.name)

    numbers = []
    taxon_sets = defaultdict(list)

    with open(taxon_set_file) as f:
        for l in f:
            number = l.strip("\n").split("\t")[0]
            taxon_set = l.strip("\n").split("\t")[1]
            
            numbers.append(number)
            
            taxon_sets[number] = taxon_set


def write_tree_model(taxon_sets, xml_file):

    for key, value in taxon_sets.items():
        name = "taxon_set_" + str(key)
        
        xml_file.write(f'<tmrcaStatistic id="tmrca({name})" absolute="false"  includeStem="false">')
        xml_file.write("\t<mrca>")
        xml_file.write(f"\t\t<taxa idref={name}/>")
        xml_file.write("\t</mrca>")
        xml_file.write('\t<treeModel idref="treeModel"/>')
        xml_file.write("</tmrcaStatistic>")
        
        xml_file.write(f'<tmrcaStatistic id="age({name})" absolute="true" includeStem="false">')
        xml_file.write("\t<mrca>")
        xml_file.write(f"\t\t<taxa idref={name}/>")
        xml_file.write("\t</mrca>")
        xml_file.write('\t<treeModel idref="treeModel"/>')
        xml_file.write("</tmrcaStatistic>")
        
        xml_file.write(f'<monophylyStatistic id="age({name})">')
        xml_file.write("\t<mrca>")
        xml_file.write(f"\t\t<taxa idref={name}/>")
        xml_file.write("\t</mrca>")
        xml_file.write('\t<treeModel idref="treeModel"/>')
        xml_file.write("</monophylyStatistic>")

def write_taxa_sets(taxon_sets, xml_file):

    for key, value in taxon_sets.items():
    
        name = "taxon_set_" + str(key)
        
        xml_file.write(f'<taxa id="{name}">')
        
        for i in value:
            xml_file.write(f'\t<taxon idref="{i}"/>')
            
        xml_file.write("</taxa>")

def write_idrefs_tree_stats(taxon_sets, xml_file):

    for i in range(1,len(taxon_sets)):
        name = "taxon_set_" + str(i)
        xml_file.write(f'<monophylyStatistic idref="monophyly({name})"/>')
    for i in range(1,len(taxon_sets)):
        name = "taxon_set_" + str(i)
        xml_file.write(f'<tmrcaStatistic idref="tmrca({name})"/>')
    for i in range(1,len(taxon_sets)):
        name = "taxon_set_" + str(i)
        xml_file.write(f'<tmcraStatistic idref="age({name})"/>')