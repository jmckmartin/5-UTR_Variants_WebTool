import gene_classes as gc
import gene_functions as gf
import sys

file = "data_files/MANE.GRCh38.v0.91.select_ensembl_rna.fna"
x = gf.read_MANE(file)
n = 25

print(x[n].gene_id)
print(x[n].transcript_id)
print(x[n].chrom)
print(x[n].gene_start_coord)
print(x[n].gene_end_coord)
print(x[n].gene_strand)
print(x[n].gene_id)
print(x[n].gene_symbol)
print(x[n].hgnc_symbol)
print(x[n].ncbi_symbol)
print(x[n].gene_description)
print(x[n].sequence)