'''this module is used to get all the header info from the OrthoFinder2 results - Single_Copy_Orthologue_Sequences
1) It takes the single-copy gene fasta as input'''
import os
from gffdb import *

def GetSGfilename(pathname):
    '''this function grep all the SingleCopy fasta filename from the OrthoFinder2 results'''
    tmp_list = []
    for gff in os.listdir(pathname):
        tmp_list.append(str(gff))
    return tmp_list


def GetsequenceID(fa):
    '''this function is used to grep the sequence id in a .fa file
    which means user have to deploy the function multiple times to get all the sequence ID in the SingleCopy dir'''
    tmp_list = []
    with open(fa, 'r') as input:
        for line in input:
            if ">"  in line:
                tmp_list.append(line.strip('>').strip('\n'))
    input.close()
    return tmp_list


def Exonbed(gffdb_dict, gffdata_path, sg_fasta_path):
    # Load the database into the environment
    bed_list = []   # build the list to deposit the bed info for each target chromosome
    gffdb_keys_list = list(gffdb_dict.keys())    # using the key to index the number of the list to deposit the bed info
    # print(gffdb_keys_list)   # show how many database we have to deal
    for i in range(len(gffdb_keys_list)):
        bed_list.append([])
    print(f"There will be {len(bed_list)} BED files")


    sc_fasta_list = GetSGfilename(sg_fasta_path)
    for fa in sc_fasta_list:
        fa_seqid_list = GetsequenceID(sg_fasta_path + '/' + fa)
        
        for seqid in fa_seqid_list:
                for x in gffdb_keys_list:
                    tmp_dbname = x.split('.')[0]
                    dbname = tmp_dbname.split('/')[1]
                    final_inputdb = gffdata_path + '/' + dbname + '.db'
                    # print(final_inputdb)
                    gffdb = Loadgffdb(final_inputdb)

                    try:
                        gene = gffdb[seqid]
                        parent_ID = gene['Parent'][0]
                        # print(type(parent_ID))
                        for i in gffdb.children(parent_ID, featuretype='exon', order_by='start'):
                            tmp_bed = i[0] + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i['ID']).strip('[]').strip("'") + '\n'
                            # print(tmp_bed)
                            bed_list[ gffdb_dict[x] ].append(tmp_bed)
                    except:
                        if dbname == 'Oryza_sativa':
                            clean_seqid = 'gene:' + seqid.split('-')[0].replace('t', 'g')
                            for i in gffdb.children(clean_seqid, featuretype='exon', order_by='start'):
                                tmp_bed = i[0] + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i['exon_id']).strip('[]').strip("'") + '\n'
                                # print(tmp_bed)
                                bed_list[ gffdb_dict[x] ].append(tmp_bed)
                            # print('The rice bed is successfully added')
                        else:
                            continue
    return bed_list

def Writebed(gffdb_dict, bed_list, path):
    gffdb_keys_list = list(gffdb_dict.keys())
    # print(gffdb_keys_list)
    for x in gffdb_keys_list:
        tmp_dbname = x.split('.')[0]
        dbname = tmp_dbname.split('/')[1]
        output_name = dbname + '.bed'
        # print(output_name)
        with open(path + '/' + output_name, 'w') as output:
            output.writelines(bed_list[ gffdb_dict[x] ])
            # print(gffdb_dict[x])
        output.close()