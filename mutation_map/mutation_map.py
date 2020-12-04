import subprocess
import re
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1_extended

GENOME_REF = 'NC_045512.2'

def generate_dummy_vcf(out_vcf_file='dummy.vcf'):
    """
    Generate artificial VCF-file with all possible SARS-CoV-2 nucleotide mutations
    """
    with open(out_vcf_file, "a") as file_:
        file_.write('##fileformat=VCFv4.2\n')
        file_.write('#CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'QUAL'+'\t'+
                    'FILTER'+'\t'+'INFO'+'\t'+'FORMAT'+'\t'+'sample\n')
        for num in range(1, 29728):
            for nuc in [('A', 'C'), ('A', 'T'), ('A', 'G'), ('C', 'A'), ('C', 'T'), ('C', 'G'),
                        ('G', 'C'), ('G', 'T'), ('G', 'A'), ('T', 'C'), ('T', 'G'), ('T', 'A')]:
                file_.write(GENOME_REF+'\t'+str(num)+'\t'+'.'+'\t'+nuc[0]+'\t'+nuc[1]+'\n')


def launch_snpeff(in_vcf_file='dummy.vcf', out_vcf_file='dummy_snpeff.vcf'):
    """
    Launch Snpeff to annotate VCF-file
    """
    cmd_parts = ['java', '-jar', 'snpEff.jar', "eff", GENOME_REF, in_vcf_file,
                 "1>", out_vcf_file]
    subprocess.run(" ".join(cmd_parts), shell=True, check=True)


def convert_protein_mutations_from_3_to_1_letters(muts: [list, set]):
    """
    Convert protein mutations from 3-letter acids to 1-letter acid format.
    Example: "p.Thr5262Ile" -> "T5262I"
    """
    new_muts = []
    for mut in muts:
        m = re.match(r"p\.(?P<acid1>[A-Z][a-z][a-z])(?P<pos>\d+)(?P<acid2>[A-Z][a-z][a-z])", mut)
        try:
            assert m, "Unexpected format (correct example: 'p.Thr42Ser')."
            acid1 = m['acid1']
            acid2 = m['acid2']
            assert acid1 in protein_letters_3to1_extended, f'Cannot recognize acid1: {acid1}'
            assert acid2 in protein_letters_3to1_extended, f'Cannot recognize acid2: {acid2}'
            new_acid1 = protein_letters_3to1_extended[acid1]
            new_acid2 = protein_letters_3to1_extended[acid2]
            new_mut = f"{new_acid1}{m['pos']}{new_acid2}"
            new_muts.append(new_mut)
        except AssertionError as e:
            #print(f"Warning while parsing mutation '{mut}' -> it will be skipped. Details: {e}")
            pass
    return new_muts

def extract_protein_mutations(info_text):
    """
    Example: QNAME=hCoV-19...;QSTART=274;QSTRAND=+;ANN=
    T|synonymous_variant|LOW|ORF1ab|GU280_gp01|transcript|GU280_gp01|
    protein_coding|1/2|c.9C>T|p.Ser3Ser|9/21291|9/21291|3/7096||,
    T|synonymous_variant|LOW|ORF1ab|GU280_gp01|transcript|YP_009725297.1|
    protein_coding|1/1|c.9C>T|p.Ser3Ser|9/540|9/540|3/179||WARNING_TRANSCRIPT_NO_STOP_CODON,
    ...,
    T|upstream_gene_variant|MODIFIER|ORF1ab|GU280_gp01|transcript|YP_009742610.1|
    protein_coding||c.-2446C>T|||||2446|WARNING_TRANSCRIPT_NO_START_CODON
    """
    # Use regexp to find nucleotide mutations (for future use) and protein mutations
    res_list = []
    for m in re.finditer(r"protein_coding\|\d+/\d+\|(?P<nuc_mut>[c.\dACGT>]*)\|(?P<prot_mut>[^|]*)", info_text):
        res_list.append(m.group('prot_mut'))
    return res_list

def make_mutation_mapping_csv(in_vcf_file_snpeff='dummy_snpeff.vcf', out_csv_file='mutation_mapping.csv'):    
    # Count number of comment lines at file beginning
    comment_lines_cnt = 0
    with open(in_vcf_file_snpeff, 'r') as f:
        for i, s in enumerate(f):
            if not s.startswith('##'):
                comment_lines_cnt = i
                break
    df_vcf = pd.read_csv(in_vcf_file_snpeff, sep='\t', skiprows=comment_lines_cnt)
    with open(out_csv_file, "a") as file_w:
        file_w.write('nucleotide' + ';' + 'protein' + '\n')
        for i, (_, row) in enumerate(df_vcf.iterrows()):
            if len(row['REF']) > 1 or row['REF'] == 'N' or len(row['ALT']) > 1 or row['ALT'] == 'N':
                # Skip non-relevant mutations
                continue
            nuc_mutation = str(row['POS']) + row['REF'] + '>' + row['ALT']

            muts = extract_protein_mutations(row['INFO'])
            muts = set(muts)
            muts = convert_protein_mutations_from_3_to_1_letters(muts)
            for mut in muts:
                file_w.write(nuc_mutation + ";" + mut + '\n')

if __name__ == '__main__':
    print('Generating dummy VCF...')
    generate_dummy_vcf()

    print('Launching Snpeff...')
    launch_snpeff()

    print('Making mutation mapping...')
    make_mutation_mapping_csv()

    print('Done!')
