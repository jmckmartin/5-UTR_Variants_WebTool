def Gene_Parse(line, tempvars): 
    #Still need to store 


    #line = ">ENST00000342066.8 cdna chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SAMD11 description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706]"
    #line.split(" ") - split line by spaces, will give\

    '''
    The following code essentially reads in the data in a header file, parses it using 
    layers of split functions, and assigns the correct info to temporary variables, which will
    be used to instantiate an Gene object
    '''
    Header_Data = line.split(" ", 7)

    Pre_Parse_Transcript_ID = Header_Data[0].split(">")
    tempvars["Temp_Transcript_ID"] = Pre_Parse_Transcript_ID[1]

    Pre_Parse_Chr_Coord_Strand = Header_Data[2].split(":")
    tempvars["Temp_Chr"] = Pre_Parse_Chr_Coord_Strand[2]
    tempvars["Temp_Gene_Start_Coord"] = Pre_Parse_Chr_Coord_Strand[3]
    tempvars["Temp_Gene_End_Coord"] = Pre_Parse_Chr_Coord_Strand[4]
    tempvars["Temp_Gene_Strand"] = Pre_Parse_Chr_Coord_Strand[5]

    Pre_Parse_Gene_ID = Header_Data[3].split(":")
    tempvars["Temp_Gene_ID"] = Pre_Parse_Gene_ID[1]

    Pre_Parse_Gene_Symbol = Header_Data[6].split(":")
    tempvars["Temp_Gene_Symbol"] = Pre_Parse_Gene_Symbol[1]

    Pre_Parse_Gene_Description, Pre_Parse_Gene_Info = Header_Data[7].split("[",1)

    Pre_Parse_HGNC_Symbol = Pre_Parse_Gene_Info.split(":")[3]
    tempvars["Temp_HGNC_Symbol"] = Pre_Parse_HGNC_Symbol[:-1] #remove end "]"

    Pre_Parse_Gene_Description = Pre_Parse_Gene_Description.split(":")[1]
    tempvars["Temp_Gene_Description"] = Pre_Parse_Gene_Description[:-1] #remove end whitespace

    return(tempvars)

def resetvars(vars_dict):
    vars_dict["Temp_Transcript_ID"] = ''
    vars_dict["Temp_Chr"] = ''
    vars_dict["Temp_Gene_Start_Coord"] = ''
    vars_dict["Temp_Gene_End_Coord"] = ''
    vars_dict["Temp_Gene_Strand"] = ''
    vars_dict["Temp_Gene_ID"] = ''
    vars_dict["Temp_Gene_Symbol"] = ''
    vars_dict["Temp_HGNC_Symbol"] = ''
    vars_dict["Temp_Gene_Description"] = ''
    vars_dict["Temp_cDNA_Sequence"] = ''
    return(vars_dict)