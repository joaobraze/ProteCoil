#! /usr/bin/env python

desc = """ProteCoil searchs for coiled coil regions in a given proteome. The input file must be a genome GenBank file. Paricoil2 tool must be installed.

"""

import os
import subprocess
import argparse
import textwrap
import re
from Bio import SeqIO
from BCBio import GFF
from collections import OrderedDict
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def fetch_info(myfile,my_p_score):

    '''Primarily, it is rised a dictionary with the coiled sequences as keys and the gene name, protein_id and product as values.
    The genome gbk file provids all the information needed. The coiled coil sequence is determined by the call of the paircoil2_printer function'''

    cc_dici={}
    cc_list=[]

    count=0
    for genome_record in SeqIO.parse(myfile, "gb"):
        for i in genome_record.features:
            if i.type == "CDS":

                for key,values in (i.qualifiers).items():
                    if key == 'translation':
                        with open("my_protein.fasta",'w') as f: #it is needed for feed the paircoil2_printer
                            f.write('>my_protein')
                            f.write('\n')
                            protein_seq=''.join(values)
                            f.write(protein_seq)
                        cc_list+=paircoil2_printer("my_protein.fasta",my_p_score) # a CC list per each protein analysed
                        count +=1
                        for cc in cc_list:
                            if cc in protein_seq:
                                try:
                                    cc_dici.setdefault(cc, []).append(''.join((i.qualifiers).get('gene'))) #some annotations don't have gene name.
                                except TypeError:
                                    cc_dici.setdefault(cc, []).append(''.join((i.qualifiers).get('locus_tag')))
                                start_pst= i.location.start + 3*protein_seq.index(cc)
                                end_pst= i.location.start + 3*(protein_seq.index(cc) + len(cc))
                                cc_dici.setdefault(cc, []).append(''.join((i.qualifiers).get('protein_id')))
                                cc_dici.setdefault(cc, []).append(''.join((i.qualifiers).get('product')))
                                cc_dici.setdefault(cc, []).append(start_pst)
                                cc_dici.setdefault(cc, []).append(end_pst)
                                cc_dici.setdefault(cc, []).append(len(protein_seq))


    '''The following lines attemp to fix appending problems related with duplicated proteins. It must be revisited in future'''


    for cc,seq_info in cc_dici.items(): #Due to duplicated proteins, certain appended info has the same info twice. The following lines eliminate the double
        if len(seq_info)==18:
            n=0
            #print (seq_info)
            for i in range(int(len(seq_info)/6)-1):
                if (seq_info[4+n]) == (seq_info[10+n]):
                    for u in range(6):
                        seq_info.pop(12+n-6)
                n+=6


    supl_dici={}    #Due to duplicated proteins, the information is appended in one CC sequence. The following lines separeted that information and mark the double sequence with an X
    supl_list=[]
    for cc,seq_info in cc_dici.items():
        if len(seq_info) == 12: #if true, means that for certain cc region duplicated information was added.
            a=seq_info.pop(6)
            b=seq_info.pop(6)
            c=seq_info.pop(6)
            e=seq_info.pop(6)
            f=seq_info.pop(6)
            g=seq_info.pop(6)
            supl_dici[(str(cc)[0:-1]+"X")]=[a,b,c,e,f,g] #the X marks the sequece as duplicated
    cc_dici.update(supl_dici)


    print ("...It was analysed %d proteins for Coiled Coil region search...\n"%count)
    print ("...Paircoil2: on this proteome %d coiled coil sequences were predicted...\n"%len(cc_list))

    return cc_dici

def paircoil2_printer(my_protein,my_p_score):

    cmd = ["paircoil2",my_protein,"paircoil2_output"]
    subprocess.run(cmd)

    cc_seq=''
    lista=[]
    with open("paircoil2_output",'r') as handle:
      for each in (handle.readlines())[4:-1]:
          try:
              p_score=float(each[37:44])
              aa=each[15]
              if p_score < my_p_score:  #set the p-score
                  cc_seq +=aa
              else:
                  if cc_seq!='':
                      lista.append(cc_seq)
                  cc_seq='' #reset the value of the variable cc_seq

          except ValueError:
              pass
              print("error in Paircoil2")

    lista.append(cc_seq) #For the rarely situations, when the predicted CC go until the end of the protein sequence

    return lista


def gff3_generator(cc_info,gbk_file):

    '''The following lines define the header of the GFF3 file'''

    record = SeqIO.read(gbk_file, "gb")
    rec = SeqRecord(record,str(gbk_file)[:-3])

    qualifiers = {"source": "RefSeq", "score": 0.0, "other": ["..."],
                  "ID":""}

    top_feature = SeqFeature(FeatureLocation(0, len(rec)), type="region", strand=1,
                             qualifiers=qualifiers)
    rec.features = [top_feature]

    '''next happens the connection to the dictionary.
    the target data is the dict.values: a list with the cc information. This list can be easily manipulated according to what it is needed'''

    nr=0
    for i in sorted(cc_info.values(), key=lambda x: x[3]):

        nr+=1
        gene_parent=i[0]
        protein_name = i[1]
        protein_product=i[2]
        startin_genome = i[3]  #start position in genome
        endin_genome = i[4]  #end position in genome


        qualifiers = {"source": "PAIRCOIL2_prediction", "score": 0.0, "ID":"j"+str(nr),
                      "Parent":gene_parent, "Protein_id":protein_name,"Product":protein_product}

        b = SeqFeature(FeatureLocation(startin_genome, endin_genome), type="CDS", strand=1,
                                 qualifiers=qualifiers)

        top_feature.sub_features.append(b)


    rec.features = [top_feature]


    '''Last step: all the previous data is written in a GFF3 file'''

    my_output=str(gbk_file)[:-4]+".gff3" #the name of my gff3 file

    with open("results/"+my_output, "w") as out_handle:
        GFF.write([rec], out_handle)

    print ("..."+ my_output + " was created")

###############################################################################
def main(gbk_file,p_score):

    cc_dici=fetch_info(gbk_file,p_score)

    gff3_generator(cc_dici,gbk_file)

    return cc_dici

##########################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("genomefile", help="the genome GenBank file name")
    parser.add_argument("-p","--p_score", type=float, default=0.1, help="the cutoff p-score")
    args = parser.parse_args()

    main(args.genomefile,args.p_score)
