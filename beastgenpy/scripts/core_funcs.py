from Bio import SeqIO
import datetime as dt

def decimal_date(date_string):

    date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    final_date = date.year + float(date.toordinal() - start) / year_length

    return str(final_date)

def write_tax_ID_blocks(fasta, traits, trait_locs, xml_file):

    names = []
    sequences = {}

    aln = SeqIO.parse(fasta, 'fasta')

    for seq in aln:
        names.append(seq.id)
        sequences[name] = seq.seq
        aln_len = len(seq.seq)

    number_seqs = len(names)

    xml_file.write(f'/t<!-- The list of taxa to be analysed (can also include dates/ages).          -->\n\t<!-- ntax={number_seqs}                                                               -->\n')
    xml_file.write('/t<taxa id="taxa">')

    for sequence in names:
        name_parts = sequence.split("|")
        dec_date = decimal_date(name_parts[-1])
        xml_file.write(f'\t\t<taxon id="{sequence}">\n')
        xml_file.write(f'\t\t\t<date={dec_date} direction="forwards" units="years"/>\n')
        for trait in traits:
            trait_value = name_parts[trait_locs[trait]]
            xml_file.write(f'\t\t\t<attr name="{trait}">\n')
            xml_file.write(f'\t\t\t\t{trait_value}\n')
            xml_file.write(f'\t\t\t</attr>\n')
        xml_file.write('\t\t</taxon>\n')

    xml_file.write("\t</taxa>\n")

    return sequences, aln_len, number_seqs
    

def write_alignment(sequences, xml_file, aln_len, number_seqs):

    xml_file.write('\t<!-- sequence alignment --> \n')
    xml_file.write(f'\t<!-- ntax={number_seqs} nchar={aln_len} -->\n')
    
    xml_file.write('\t<alignment id="alignment" dataType="nucleotide">\n')


    for name, seq in sequences.items():
        xml_file.write('\t\t<sequence>\n')
        xml_file.write(f'\t\t\t<taxon idref="{name}"/>\n')
        xml_file.write(f'\t\t\t{seq}\n')
        xml_file.write('\t\t</sequence>\n')

    xml_file.write('\t</alignment>')

def codon_positions(xml_file, codon_partitioning):

    if not codon_partitioning:
        xml_file.write("\t<!-- Unique patterns from 1 to end -->\n")
        xml_file.write('\t<patterns id="patterns" from="1" strip="false"\n')
        xml_file.write('\t\t<alignment idref="alignment"/>\n')
        xml_file.write('\t</patterns>')

    ## add in here blocks for doing codon partitioning
    else:
        pass





        

    




