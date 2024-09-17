import argparse
import csv
from typing import Dict, Set


parser = argparse.ArgumentParser(description='Convert GOHierarchy to uniprot format')
parser.add_argument('input_file', type=str, help='Input file')
parser.add_argument('fasta_file', type=str, help='FASTA file to source sequence')
parser.add_argument('organism_dat', type=str, help='organism.dat file to source taxon ID')
parser.add_argument('-m', '--id_forward_mapping', type=str, help='ID forward mapping file, e.g., librarySeqMap')
parser.add_argument('-i', '--interpro_file', type=str, help='InterPro TSV file to source InterPro annotations')


EXP_CODES = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"]
EXCLUDED_TERMS = ["GO:0005488", "GO:0005515", "GO:0003674", "GO:0008150", "GO:0005575"]


def parse_organism_dat(organism_dat: str) -> Dict[str, str]:
    org_tax_lkp = {}
    with open(organism_dat, 'r') as org_f:
        org_reader = csv.reader(org_f, delimiter="\t")
        taxon_id = None
        for row in org_reader:
            os_code = row[2]
            taxon_id = row[5]
            if os_code:
                org_tax_lkp[os_code] = taxon_id
    return org_tax_lkp


def parse_fasta_file(fasta_file: str) -> Dict[str, str]:
    fasta_seq_lkp = {}
    seq_id = ''
    with open(fasta_file, 'r') as fasta_f:
        for l in fasta_f:
            if l.startswith('>'):
                seq_id = l.split('>', maxsplit=1)[1].strip()
                fasta_seq_lkp[seq_id] = ''
            else:
                fasta_seq_lkp[seq_id] += l.strip()
    return fasta_seq_lkp


def parse_id_forward_mapping(id_forward_mapping_file: str) -> Dict[str, str]:
    id_fwd_lkp = {}
    with open(id_forward_mapping_file, 'r') as id_fwd_f:
        reader = csv.reader(id_fwd_f, delimiter='\t')
        for row in reader:
            prev_id = row[0]
            new_id = row[1]
            id_fwd_lkp[prev_id] = new_id
    return id_fwd_lkp


def parse_interpro_file(interpro_file: str) -> Dict[str, Set]:
    interpro_lkp = {}
    with open(interpro_file, 'r') as interpro_f:
        reader = csv.reader(interpro_f, delimiter='\t')
        for r in reader:
            uniprot_id = r[0]
            interpro_id = r[1]
            if uniprot_id not in interpro_lkp:
                interpro_lkp[uniprot_id] = set()
            interpro_lkp[uniprot_id].add(interpro_id)
    return interpro_lkp


if __name__ == "__main__":
    args = parser.parse_args()

    organism_taxon_lkp = parse_organism_dat(args.organism_dat)
    fasta_sequence_lkp = parse_fasta_file(args.fasta_file)
    id_forward_mapping = None
    if args.id_forward_mapping:
        id_forward_mapping = parse_id_forward_mapping(args.id_forward_mapping)
    interpros = {}
    if args.interpro_file:
        interpros = parse_interpro_file(args.interpro_file)

    with open(args.input_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')

        genes = {}
        for row in reader:
            sequence_id = row['SequenceID']
            go_term = row['GOHierarchy'].split('>')[0]
            evidence_code = row['EvidenceCode']
            with_info = row['With'] if row['With'] else ''
            reference = row['Reference']
            date = row['Date']
            db = row['DB']

            if evidence_code not in EXP_CODES:
                continue

            if go_term in EXCLUDED_TERMS:
                continue

            if id_forward_mapping:
                sequence_id = id_forward_mapping.get(sequence_id)
                if sequence_id is None:
                    # Skip if sequence_id is not found in the id_forward_mapping
                    continue

            if sequence_id not in genes:
                genes[sequence_id] = set()
            genes[sequence_id].add((go_term, evidence_code))

        for sequence_id, go_annots in genes.items():
            os_code, gene, uniprot = sequence_id.split('|')
            uniprot_id = uniprot.split('=')[-1]
            taxon_id = organism_taxon_lkp[os_code]
            print(f"ID   {sequence_id}")
            print(f"AC   {uniprot_id};")
            print(f"OX   NCBI_TaxID={taxon_id};")
            for go_term, evidence_code in go_annots:
                # Ex: DR   GO; GO:0022625; C:cytosolic large ribosomal subunit; IEA:Ensembl.
                print(f"DR   GO; {go_term}; C:fake term label; {evidence_code}:FakeDB.")
            if uniprot_id in interpros:
                for interpro_id in interpros[uniprot_id]:
                    print(f"DR   InterPro; {interpro_id}; Fake_interpro_label.")
            sequence = fasta_sequence_lkp.get(sequence_id)
            if sequence:
                print(f"SQ   SEQUENCE header placeholder")
                print(f"SQ   {sequence}")
            print(f"//")
