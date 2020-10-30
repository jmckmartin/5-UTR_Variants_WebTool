def read_MANE(filename):
    Sequence = ''
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>'
            '''
            if temporary values not empty, assign using:
                Create object (with):
                    temporary_gene_name = Gene(self, temporary values)
            
            Reset Temporary variables
            Sequence = ''
            Parse data and assign into temporary variables
            '''
            if not line[0] == '>':
                Sequence += line.rstrip()

#>ENST00000342066.8 cdna chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SAMD11 description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706]

def Gene_Parse(): #Cumbersome version of Parsing, further break into functions, rather than re-using code
line.split(" ") - split line by spaces, will give\

line = ">ENST00000342066.8 cdna chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SAMD11 description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706]"

'''
The following code essentially reads in the data in a header file, parses it using 
layers of split functions, and assigns the correct info to temporary variables, which will
be used to instantiate an Gene object
'''

Header_Data = line.split(" ", 7)

Pre_Parse_Transcript_ID = Header_Data[0].split(">")
Temp_Transcript_ID = Pre_Parse_Transcript_ID[1]

Pre_Parse_Chr_Coord_Strand = Header_Data[2].split(":")
Temp_Chr = Pre_Parse_Chr_Coord_Strand[2]
Temp_Gene_Start_Coord = Pre_Parse_Chr_Coord_Strand[3]
Temp_Gene_End_Coord = Pre_Parse_Chr_Coord_Strand[4]
Temp_Gene_Strand = Pre_Parse_Chr_Coord_Strand[5]

Pre_Parse_Gene_ID = Header_Data[3].split(":")
Temp_Gene_ID = Pre_Parse_Gene_ID[1]

Pre_Parse_Gene_Symbol = Header_Data[6].split(":")
Temp_Gene_Symbol = Pre_Parse_Gene_Symbol[1]

Pre_Parse_Gene_Description, Pre_Parse_Gene_Info = Header_Data[7].split("[",1)

Pre_Parse_HGNC_Symbol = Pre_Parse_Gene_Info.split(":")[3]
Temp_HGNC_Symbol = Pre_Parse_HGNC_Symbol[:-1] #remove end "]"

Pre_Parse_Gene_Description = Pre_Parse_Gene_Description.split(":")[1]
Temp_Gene_Description = Pre_Parse_Gene_Description[:-1] #remove end whitespace



class Gene:
    def __init__(self, ):
        self.Transcript_ID = 
        self.Chromosome = 
        self.Genetic_Coordinates = 
        self.Gene_ID = 
        self.Symbol = 
        self.Description =
        self.Hgnc_ID = 
        self.cDNA_Sequence = 

    '''
    Other Functions :: Get 5'UTR Sequence
    '''
    