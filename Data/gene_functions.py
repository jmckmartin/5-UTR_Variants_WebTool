from gene_classes import Gene

def reset_vars(vars_dict):
    vars_dict["Temp_Transcript_ID"] = ''
    vars_dict["Temp_Chr"] = ''
    vars_dict["Temp_Gene_Start_Coord"] = ''
    vars_dict["Temp_Gene_End_Coord"] = ''
    vars_dict["Temp_Gene_Strand"] = ''
    vars_dict["Temp_Gene_ID"] = ''
    vars_dict["Temp_Gene_Symbol"] = ''
    vars_dict["Temp_HGNC_Symbol"] = ''
    vars_dict["Temp_NCBI_Symbol"] = ''
    vars_dict["Temp_Gene_Description"] = ''
    vars_dict["Temp_cDNA_Sequence"] = ''
    return vars_dict


def read_MANE(filename):
    FirstEntry = True

    # Temporary variables
    tempvars = {
        "Temp_Transcript_ID": '',
        "Temp_Chr": '',
        "Temp_Gene_Start_Coord": '',
        "Temp_Gene_End_Coord": '',
        "Temp_Gene_Strand": '',
        "Temp_Gene_ID": '',
        "Temp_Gene_Symbol": '',
        "Temp_HGNC_Symbol": '',
        "Temp_NCBI_Symbol": '',
        "Temp_Gene_Description": '',
        "Temp_cDNA_Sequence": ''
    }

    Genelist = []
    # Stores all Gene objects as the file is read

    with open(filename, 'r') as f:
        # for-else clause executed when the loop terminates through exhaustion of the list
        gene_count = 0
        for line in f:
            if line[0] == '>':  # ie, new gene reached

                if (FirstEntry == True):
                    # If this is the first time a gene is read, then no previous information collected, Parse header and proceed

                    gene_parse(line, tempvars)
                    FirstEntry = False  # Only reached after first gene recognised, this in future iterations will pass through if statement 2

                else: #elif(FirstEntry == False):  # Note, elif here?
                    # All data, including sequence for previous data stored in temp variables, thus instantiate object
                    # with this data, reset, and begin to store information for new gene
                    gene_count += 1
                    print(gene_count)
                    print(tempvars["Temp_Gene_Symbol"])
                    tempvars["Temp_Gene_Symbol"] = Gene(tempvars)
                    # Construct Object, initialising with temporary variables
                    Genelist.append(tempvars["Temp_Gene_Symbol"])
                    # - Add object to data structure

                    reset_vars(tempvars)  # - Reset Variables
                    gene_parse(line, tempvars)  # - Parse Header and assign new variables


            if line[0] != '>':
                tempvars["Temp_cDNA_Sequence"] += line.rstrip()

        else:
            gene_count += 1
            print(gene_count)
            print(tempvars["Temp_Gene_Symbol"])
            tempvars["Temp_Gene_Symbol"] = Gene(tempvars)
            Genelist.append(tempvars["Temp_Gene_Symbol"])
            reset_vars(tempvars)

    return(Genelist)

def gene_parse(line, tempvars):

    #Add formal definition so can be found in "help"
    '''
    The following code essentially reads in the data in a header file, parses it using
    layers of split functions, and assigns the correct info to temporary variables, which will
    be used to instantiate an Gene object
    '''

    Header_Data = line.split(" ", 7)
    # Header Data - Splits Header into 8 chunks, where the splits occur at ' '

    Pre_Parse_Transcript_ID = Header_Data[0].split(">")
    # Pre_Parse_Transript_ID - Splits first chunk (eg >ENST00000342066.8) into two, sepeperated by '>', and takes the 2nd element
    tempvars["Temp_Transcript_ID"] = Pre_Parse_Transcript_ID[1]

    Pre_Parse_Chr_Coord_Strand = Header_Data[2].split(":")
    # Pre_Parse_Chr_Coord_Strand - Further splits third "Header Data" chunk (chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12)
    # into smaller parts, seperated by ":" - these parts are then assigned to the relevant information
    tempvars["Temp_Chr"] = Pre_Parse_Chr_Coord_Strand[2]
    tempvars["Temp_Gene_Start_Coord"] = Pre_Parse_Chr_Coord_Strand[3]
    tempvars["Temp_Gene_End_Coord"] = Pre_Parse_Chr_Coord_Strand[4]
    tempvars["Temp_Gene_Strand"] = Pre_Parse_Chr_Coord_Strand[5]

    Pre_Parse_Gene_ID = Header_Data[3].split(":")
    # Pre_Parse_Gene_ID - Further splits fourth "Header Data" chunk (eg gene:ENSG00000187634.12) into two,
    # seperated by ":", where the second element is used for the gene ID
    tempvars["Temp_Gene_ID"] = Pre_Parse_Gene_ID[1]

    Pre_Parse_Gene_Symbol = Header_Data[6].split(":")
    # Pre_Parse_Gene_Symbol = splits seventh chunk (eg gene_symbol:SAMD11), across ":", and so symbol can be assigned
    tempvars["Temp_Gene_Symbol"] = Pre_Parse_Gene_Symbol[1]
    # if Pre_Parse_Gene_Symbol[1] == "SMIM40":
    #     print ("TEST NUMBER")
    #     sys.exit()


    Pre_Parse_Gene_Description, Pre_Parse_Gene_Info = Header_Data[7].split("[", 1)
    # Pre_Parse_Gene_Description = splits eighth chunk (eg description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706])
    # into 2, seperated by "[", takes the first element, further splits by ":", and takes the second element, before cutting off whitespace

    # Pre_Parse_Gene_Info = splits eighth chunk (eg description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706])
    # into 2, seperated by "[", takes the second element, further splits according to ":", takes the last element, and removes the end "]"

    Pre_Parse_Symbol_Check = Pre_Parse_Gene_Info.split(" ", 1)[0]
    # Pre_Parse_Symbol_Check - splits up chunk in order to see the source of the Gene information as HGNC / NCBI

    if Pre_Parse_Symbol_Check == "Source:HGNC":
        # If the source is HGNC, this block is entered, the test is split, and temporary variables are assigned
        Pre_Parse_HGNC_Symbol = Pre_Parse_Gene_Info.split(":")[3]
        tempvars["Temp_HGNC_Symbol"] = "HGNC:" + Pre_Parse_HGNC_Symbol[:-2]  # remove end "]"

    elif Pre_Parse_Symbol_Check == "Source:NCBI":
        # If the source is NCBI, this block is entered, the test is split, and temporary variables are assigned
        Pre_Parse_NCBI_Symbol = Pre_Parse_Gene_Info.split(":")[2]
        tempvars["Temp_NCBI_Symbol"] = Pre_Parse_NCBI_Symbol[:-2] # remove end "]"

    Pre_Parse_Gene_Description = Pre_Parse_Gene_Description.split(":")[1]
    tempvars["Temp_Gene_Description"] = Pre_Parse_Gene_Description[:-1]  # remove end whitespace

    return(tempvars)