""" 
The form of the data in the "MANE.GRCh38.v0.91.select_ensembl_rna.fna" file is 
>ENST00000342066.8 cdna chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SAMD11 description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706]
$$sequence$$
>
Structure is thus:
> | Transcript ID (ENST...) | cdna | chromosome_label:Assembly:start:end:strand | gene_label:ID (ENSG) | genebiotype_label:type | transcript_biotype:type | symbol_label:symbol | description_label:string(incl spaces) | [Source_label:source | Symbol_label...]

-It seems then that the best way would be to recognise new Transcript by >
- Relevant information could then be sorted using count for spaces and colons 
"""