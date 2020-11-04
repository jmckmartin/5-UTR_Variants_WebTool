def read_MANE(filename):
    
    FirstEntry = True

    #Temporary variables
    tempvars = {
    "Temp_Transcript_ID" : '',
    "Temp_Chr" : '',
    "Temp_Gene_Start_Coord" : '',
    "Temp_Gene_End_Coord" : '',
    "Temp_Gene_Strand" : '',
    "Temp_Gene_ID" : '',
    "Temp_Gene_Symbol" : '',
    "Temp_HGNC_Symbol" : '',
    "Temp_Gene_Description" : '',
    "Temp_cDNA_Sequence" : ''
    }

    with open("Data/MANE.GRCh38.v0.91.select_ensembl_rna.fna", 'r') as f:
        #for-else clause executed when the loop terminates through exhaustion of the list
        for line in f:

            if line[0] == '>': #ie, new gene reached
                
                if FirstEntry == True: 
                    #If this is the first time a gene is read, then no previous information collected, Parse header and proceed
                    
                    #Gene_Parse(line, tempvars)

                if FirstEntry == False: #Note, elif here?
                    #All data, including sequence for previous data stored in temp variables, thus instantiate object
                    #with this data, reset, and begin to store information for new gene

                    #Construct Object, initialising with temporary variables

                    #- Add object to data structure

                    #resetvars(tempvars)# - Reset Variables

                    #Gene_Parse(line, tempvars) #- Parse Header and assign new variables

                FirstEntry == False #Only reached after first gene recognised, this in future iterations will pass through if statement 2
                
            if line[0] != '>':
                tempvars["Temp_cDNA_Sequence"] += line.rstrip()

        else: 
            #Construct Object, initialising with temporary variables
            #- Add object to data structure

            #resetvars(tempvars)# - Reset Variables
