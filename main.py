'''this is main workspace of Jtools'''
from gffdb import *
from singlecopyFinder import *
from multiprocessing import Process
import sys

'''Using sys to parse the gff input in the commandline mode'''
# gffname_list = []
# def ParseGFF(*args):
#     for index, name in enumerate(sys.argv):
#         if index != 0:
#             gffname_list.append(name)
# ParseGFF()
# print(gffname_list)


'''test'''



'''Run the Makegffdb function on the multiprocessiing mode'''
if __name__ == '__main__':
    '''Note: the variables'''
    gffpath = 'gffpath'
    gffdata_path = 'gffdatabase'
    sg_fasta_path = 'Single_Copy_Orthologue_Sequences'
    bedoutput_patn = 'BED_of_Exon_of_SingleCopygene'

    # Make GFFdatabase according to the specified GFF name list
    def MultirunMakeGFFdb(tmp_list):
        '''tmp_list is a list which contains the input gff name'''
        process_list = []
        for gff in tmp_list:
            p = Process(target=Makegffdb, args=(gff,))
            p.start()
            process_list.append(p)
        for p in process_list:
            p.join()
    gffname_list = Readgff(gffpath)
    print(gffname_list)
    # gffname_list = ['moso.hic.chr.gff', 'sub_A1_coding_gene.gff3', 'sub_A2_coding_gene.gff3', 'sub_B1_coding_gene.gff3', 'sub_B2_coding_gene.gff3', 'sub_C1_coding_gene.gff3', 'sub_C2_coding_gene.gff3', 'Oryza_sativa.IRGSP-1.0.44.gff3.gz']
    # the gffname_list could be "argparse" in put later
    order_list = [ x for x in range(0, len(gffname_list) + 1) if x!=8]  # make sure the init number is 1 (1-based order) 
    # print(order_list)
    gffdb_dict = dict(zip(gffname_list, order_list))
    # print(gffdb_dict)
    MultirunMakeGFFdb(gffname_list)
    

    # BED file generation
    bed_list = Exonbed(gffdb_dict, gffdata_path, sg_fasta_path)
    Writebed(gffdb_dict, bed_list, bedoutput_patn)

    '''Scripts'''
    # Now We hava gffutils feature database: moso, sub_A1_coding_gene, sub_A2_coding_gene, sub_B1_coding_gene, sub_B2_coding_gene, sub_C1_coding_gene, sub_C1_coding_gene
    # Ped analysis
    # mosodb = Loadgffdb('gffdatabase/moso.db')
    # print(type(mosodb))
    # sc_fasta_list = GetSGfilename(sg_fasta_path)
    # # Ped14_exon_list = []
    # Ped14_exon_bed = []
    # # Ped16_exon_list = []
    # Ped16_exon_bed = []
    # for fa in sc_fasta_list:
    #     tmp_seqid_list = GetsequenceID(sg_fasta_path + '/' + fa)
    #     for index, seqid in enumerate(tmp_seqid_list):
    #         # parse_seqid = seqid.split('.')[0]
    #         try:
    #             mosodb[seqid]
    #         except:
    #             continue
    #         if index == 0:
    #             for i in mosodb.children(seqid, featuretype='exon', order_by='start'):
    #                 tmp_bed = i[0] + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i['ID']).strip('[]').strip("'") + '\n'
    #                 Ped14_exon_bed.append(tmp_bed)

    #         elif index == 1:
    #             for i in mosodb.children(seqid, featuretype='exon', order_by='start'):
    #                 tmp_bed = i[0] + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i['ID']).strip('[]').strip("'") + '\n'
    #                 Ped16_exon_bed.append(tmp_bed)
    # # print(len((Ped14_exon_list)))
    # # print(len(Ped16_exon_list))
    # # print(Ped14_exon_list[0])
    # # print(Ped14_exon_bed[0])

    # # Write Exon BED files
    # with open('BED_of_Exon_of_SingleCopygene/Ped14_exon.bed', 'w') as output:
    #     output.writelines(Ped14_exon_bed)
    # output.close()

    # with open('BED_of_Exon_of_SingleCopygene/Ped16_exon.bed', 'w') as output:
    #     output.writelines(Ped16_exon_bed)
    # output.close()