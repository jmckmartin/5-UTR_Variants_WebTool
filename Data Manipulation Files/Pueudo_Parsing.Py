""" 
The form of the data in the "MANE.GRCh38.v0.91.select_ensembl_rna.fna" file is 
>ENST00000342066.8 cdna chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SAMD11 description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706]
$$sequence$$
>
Structure is thus:
> | Transcript ID (ENST...) | cdna | chromosome_label:Assembly:start:end:strand | gene_label:ID (ENSG) | genebiotype_label:type | transcript_biotype:type | symbol_label:symbol | description_label:string(incl spaces) | [Source_label:source | Symbol_label...]

-It seems then that the best way would be to recognise new Transcript by >
- Relevant information could then be sorted using count for spaces and colons 

---

Rough plan:

Read file line by line
- store the lines that don't begin with an arrow in a file for sequence

-Create Gene Class Definition
-Create Parsing and Assignment Function 
---
1 - Import Class and Function Definitions
2 - Open File 
3 - Read in file
4 - When 1st character in line = 0:
- trigger parsing function to store data in header line to temporary variables
5 - For subsequent lines until you reach the next ">", read the file (minus whitespace etc) as a string to sequence variable
6 - When you hit the next '>', create class, assign temporary variables and restart until the end of the file


Note on storage - Potentially can use 


File in

Build data structure

Class definition for Gene

if line begins with '>', 

create object for Gene, containing information for:
        - transcript ID
        - Chromostome
        - Genetic Co-ordinates
        - Gene ID
        - Symbol
        - Description
        - HGNC ID
        -*Sequence*, see below

    for subsequent lines, until you reach the next ">", read sequence into data structure

    Store the created object in a data structure containing all the Genes
"""
