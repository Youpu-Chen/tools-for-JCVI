'''this module is based on gffutils to make the local gff datebase'''
import gffutils
import os 


def Readgff(gffpath):
    '''this function is used to grep the original gff file name
    the element inside this list is like below
    gffpath/moso.hic.chr.gff '''
    tmp_list = []
    for gff in os.listdir(gffpath):
        tmp_list.append('gffpath' + '/' + str(gff))
    return tmp_list


def Makegffdb(gff):
    fn = gff
    tmp_dbname = fn.split('.')[0] + '.db'   # use '.' to parse the gff file and take the first segment of character as the name of the database
    dbname = tmp_dbname.split('/')[1]
    # print(dbname)
    # path_dbname = pathname + '/' + dbname
    path_dbname = 'gffdatabase' + '/' + dbname
    # print(path_dbname)
    if os.path.exists(path_dbname) == False:     # to make sure not to overwrite the existing gffdb, because it takes a long time to build the database
        return gffutils.create_db(fn, dbfn=path_dbname, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

def Loadgffdb(gffdb):
    return gffutils.FeatureDB(gffdb, keep_order=True)

'''Scripts
This part is waiting to be developed, e.g. check the bed is empty or not'''
    # # Check whether the BED file is empty
    # check_bed = []
    # tmp_all = os.listdir('./')
    # for i in tmp_all:
    #     if os.path.splitext(i)[0] == '.bed':
    #         check_bed.append(os.path.splitext(i)[0])