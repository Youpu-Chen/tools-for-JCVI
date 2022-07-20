'''this modue is used to extract exon-to-exon structure from the blast results'''
import pandas as pd
import os
import numpy as np
import gffutils
from collections import defaultdict


def RBHexon2exon(blastA, blastB):
    '''blast1 is the A_vs_B blast file, and the blasdt2 is the B_vs_A_blast file
    Note: 
    1) the input blast file should be filtered by the JCVI "jcvi.formats.blast best" module with "best 1" parameter 
    - the input is the absolute path, thus should be parsed carefully
    2) this function could be revised to generate the raw RBH ID datasets.'''
    nameA = blastA.split('.')[0].split('/')[1]
    blastA_maxhsp1_maxts5 = pd.read_table(blastA, index_col=False, header=None)
    blastA_maxhsp1_maxts5_ID = blastA_maxhsp1_maxts5.iloc[:, [0,1]]
    blastA_maxhsp1_maxts5_ID.columns = ['A', 'B']

    nameB = blastB.split('.')[0].split('/')[1]
    blastB_maxhsp1_maxts5 = pd.read_table(blastB, index_col=False, header=None)
    blastB_maxhsp1_maxts5_ID = blastB_maxhsp1_maxts5.iloc[:, [0,1]]
    blastB_maxhsp1_maxts5_ID.columns = ['B', 'A']

    # Merge the blast results
    AB_maxhsp1_maxts5_ID = pd.merge(blastA_maxhsp1_maxts5_ID, blastB_maxhsp1_maxts5_ID, left_on="B", right_on="B", 
                                    how = "inner", suffixes=('_A', '_B'))
    # Generate the judgements
    jud = AB_maxhsp1_maxts5_ID.iloc[:,0] == AB_maxhsp1_maxts5_ID.iloc[:,2]
    RBH_AB_maxhsp1_maxts5_ID = AB_maxhsp1_maxts5_ID[jud]
    RBH_AB_maxhsp1_maxts5_ID = RBH_AB_maxhsp1_maxts5_ID.iloc[:,[0,1]] # if using the [jud] and column indexes together (e.g. [jud].iloc[:,[0,1]]), replicates will not be removed.
    RBH_AB_maxhsp1_maxts5_ID.to_csv(f'{nameA}.{nameB}.RBHID.txt', sep='\t', header=False, index=False)
    # print(RBH_AB_maxhsp1_maxts5_ID.columns)

    # If saved the df to the workingspace, then using the combinatino of ID to read them,
    # it might be easier to handle for the function StructureJudge
    return RBH_AB_maxhsp1_maxts5_ID



def StructureJudge(nameA, nameB):
    '''this function is used to check whether the RBH exons has the exon-intron-exon structure,
    Note: the output file dir in version is exon2exon
    
    There are several condtions to be met to get the exon-intron-exon
    1) the freq of exon parent (e.g. transcript) should > 1
    2) the exon ID should be sorted according to their number
    3) the corresponding gap should be the same 
    The √ form
    - exon 1 of gene 1 on ch1 of speciesA -- exon 1 of gene 1 on ch1 of speciesB
    - exon 2 of gene 1 on ch1 of speciesA -- exon 2 of gene 1 on ch1 of speciesB
    The × form (should be filted)
    - exon 1 of gene 1 on ch1 of speciesA -- exon 1 of gene 1 on ch1 of speciesB
    - exon 2 of gene 1 on ch1 of speciesA -- exon 3 of gene 1 on ch1 of speciesB (the correspoding gap should be the same)
'''

    prefix = f'{nameA}.{nameB}'
    file = prefix + '.RBHID.txt'
    if 'Ped' in nameA:
        gffdatabase1 = gffutils.FeatureDB(f'gffdatabase/Ped.db', keep_order=True)
    else:
        gffdatabase1 = gffutils.FeatureDB(f'gffdatabase/{nameA}.db', keep_order=True)
    
    gffdatabase2 = gffutils.FeatureDB(f'gffdatabase/{nameB}.db', keep_order=True)

    # Create list and dict to deposit the info
    exon1_list = []
    exon2_list = []
    exon1_freq_dict = {}
    exon2_freq_dict = {}
    freqfilter_exon1_list = []
    freqfilter_exon2_list = []
    
    exon1_transcript1_dict = defaultdict(list)
    exon12exon2_dict = {}
    exon22exon1_dict = {}


    # Filter preparation
    # file is the {species_nameA}.{species_nameB}.RBHID.txt file
    with open(file, 'r') as input:
        for line in input:
            parsed_line = line.strip('\n').split('\t')
            # saved variables
            exon1_id = parsed_line[0]
            exon2_id = parsed_line[1]
            exon1_list.append(exon1_id)
            exon2_list.append(exon2_id)
            exon12exon2_dict[exon1_id] = exon2_id
            exon22exon1_dict[exon2_id] = exon1_id
    # print(exon22exon1_dict)
    print(f"The original exon1 list have {len(exon1_list)} exons")
    input.close()

    # Count the frequency of exon-corresponding genes
    for x in exon1_list:
        exon1_line = gffdatabase1[x]
        transcript1_id = exon1_line['Parent'][0]
        exon1_freq_dict[transcript1_id] = exon1_freq_dict.get(transcript1_id, 0) + 1
    
    for x in exon2_list:
        clean_seqid = 'transcript:' + x.split('-')[0] + '-' + x.split('-')[1]
        # print(clean_seqid)
        exon2_line = gffdatabase2[clean_seqid]
        transcript2_id = exon2_line['ID'][0]
        exon2_freq_dict[transcript2_id] = exon2_freq_dict.get(transcript2_id, 0) + 1

    # 1.1) Filter the exon1 whose freq is equal to 1
    for x in exon1_list:
        # exon_id = x.strip('\n').split('\t')[0]
        exon1_line = gffdatabase1[x]
        tmp_transcript_id = exon1_line['Parent'][0]
        tmp_count = exon1_freq_dict[tmp_transcript_id]
        # print(tmp_count)
        if tmp_count == 1:
            continue
        else:
            freqfilter_exon1_list.append(x)
    
    # 1.2) Filter the exon2 whose freq is equal to 2
    for x in exon2_list:
        clean_seqid = 'transcript:' + x.split('-')[0] + '-' + x.split('-')[1]
        # print(clean_seqid)
        exon2_line = gffdatabase2[clean_seqid]
        # print(exon2_line)
        tmp_transcript_id = exon2_line['ID'][0]
        tmp_count = exon2_freq_dict[tmp_transcript_id]
        # print(tmp_count)
        if tmp_count == 1:
            continue
        else:
            freqfilter_exon2_list.append(x)
    # print(len(freqfilter_exon2_list))
    # print(freqfilter_exon1_list)
    # print(f"The freqfilter exon1 list have {len(freqfilter_exon1_list)} exons")
    # print(f'After frequency filtering, there are {len(freqfilter_exon1_list)} exons in the dataset')

    # Using the centor species exon list to filter the exon list of the species used to blast
    exon22exon1_dict_values = exon22exon1_dict.values()
    freqfilter_exon1_list = list(set(freqfilter_exon1_list).intersection(exon22exon1_dict_values))
    freqfilter_exon1_list = sorted(freqfilter_exon1_list)
    print(freqfilter_exon1_list)
    print(f"The freqfilter exon1 list have {len(freqfilter_exon1_list)} exons")

    # 2) Reorder the exon using their ID
    for x in freqfilter_exon1_list:
        # in this case the exon id is like: transcript_id + "exon" + order number
        # Note: when applied in different dataset, it should be changed.
        # print(x)
        transcript1_id = x.split('exon')[0]
        # print(transcript1_id)
        order = x.split('exon')[1]
        exon1_transcript1_dict[transcript1_id].append(int(order))
    # print(exon1_transcript1_dict)
    exon1_transcript1_dict_keys = exon1_transcript1_dict.keys()
    # Sort the exon ID
    for x in exon1_transcript1_dict_keys:
        exon1_transcript1_dict[x] = sorted(exon1_transcript1_dict[x])

    # 3) Filter the exon-exon with no same gap
    # Waiting to be developed

    # Generate the filtered RBH.txt
    with open('exon2exon/filter' + '.' + file, 'w') as output:
        for i in exon1_transcript1_dict_keys:
            exon1_count = len(exon1_transcript1_dict[i])
            for x in range(exon1_count):
                exon1_id = i + 'exon' + str(exon1_transcript1_dict[i][x])
                exon2_id = exon12exon2_dict[exon1_id]
                output.write(f'{exon1_id}\t{exon2_id}\n')
    output.close()
    print(f'The reorder is finished!')


def Getbest(path):
    '''this function is used to generate a list which contains the best blast output (generated by JCVI) '''
    bestblast_list = os.listdir(path)
    return bestblast_list


def Combination(bestblast_list):
    '''this fuction is used to generate a dict which has the key representing the species A, and the value representing the species B'''
    prefix_dict = defaultdict(list)
    for x in bestblast_list:
        seperator = '.'
        prefix1 = x.split(seperator)[0]
        prefix2 = x.split(seperator)[1]
        prefix_dict[prefix1].append(prefix2)
    return prefix_dict


def LoopRun(prefix_dict):
    for x in prefix_dict.keys():
        # The centor of RBH (main species used to blast against, in this case, the rice) should be left out.
        if len(prefix_dict[x]) >= 2:
            continue
        else:
            best1 = blastpath + '/' + x + '.' + prefix_dict[x][0] + '.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best'
            best2 = blastpath + '/' + prefix_dict[x][0] + '.' + x + '.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best'
            # print(best1)
            # print(best2)
            RBHexon2exon(best1, best2)
            StructureJudge(x, prefix_dict[x][0])


if __name__ == "__main__":
    '''Script'''
    # RBHexon2exon('Ped14.osa.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best', 'osa.Ped14.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best')
    # RBHexon2exon('Ped14.osa.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best', 'osa.Ped14.omt6.maxhsp1.maxts5.E5.id70.blast.txt.best')
    # # print(raw_RBH_df['A_A'].str.contains('\.'))
    # # print(raw_RBH_df['B'].str.contains('\.'))

    # StructureJudge('Ped14', 'osa')

    '''Serious'''
    blastpath = 'Ped14_Ped16_best'
    bestblast_list = Getbest(blastpath)
    # # print(bestblast_list)
    prefix_dict = Combination(bestblast_list)
    # # print(prefix_dict)
    # print(prefix_dict.keys())

    # # The loop below is used to run the RBHexon2exon and StructureJudge function
    LoopRun(prefix_dict)

