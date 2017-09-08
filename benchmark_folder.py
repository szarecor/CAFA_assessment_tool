#!/home/nzhou/anaconda2/python
# -*- coding: utf-8 -*-

import ID_conversion.ID_conversion
import os
import argparse
import gzip

def taxon_name_converter(taxonID):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI',
    '7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE',
    '237561':'CANAX','559292':'YEAST','284812':'SCHPO','8355':'XENLA','224308':'BACSU',
    '99287':'SALTY','243232':'METJA','321314':'SALCH','160488':'PSEPK','223283':'PSESM',
    '85962':'HELPY','243273':'MYCGE','170187':'STRPN','273057':'SULSO'}
    return taxonTable[taxonID]    
      
#The following species have their own ID mapping file in uniprot 
file_mapping_available = ['3702',
 '6239',
 '9031',
 '7955',
 '44689',
 '7227',
 '83333',
 '9606',
 '10090',
 '10116',
 '284812',
 '559292']

#For the rest of the species, we use the uniprot web conversion tool

def write_uniprot_idMapping(idmappingfile,taxon):
    #This writes the lines in uniprot id mapping file that maps accession to ID
    #File written still too large to store in memory
    #NOT IN USE
    outhandle = open('./uniprot_ac_to_id_'+str(taxon)+'.map','w')
    for line in gzip(idmappingfile,'rb'):
        fields = line.strip().split('\t')
        if fields[1]=='UniProtKB-ID':
            outhandle.write(line)
    outhandle.close()
    print('Done!\n')
    
def read_uniprot_IdMapping(idmappingfile, taxon):
    #Run only once please
    #This function returns a dictionary
    #Run species specific file please!
    if taxon in file_mapping_available:
        uniprotdict = dict()
        handle = open(idmappingfile,'rb')
        for line in handle:
            fields = line.strip().split('\t')
            uniprotdict[fields[0]]=fields[2]
                #print(line)
                #break
        handle.close()
        print('Uniprot Id Mapping read for %s.\n' % taxon )
        
    else:
        print('Species-specific file unavailable!\n')
        print('Run web conversion tool first\n')
        #uniprotdict=None
    return(uniprotdict)
    
def read_uniprot_IdMapping_fromWeb(idmappingfile,taxon):
    uniprotdict = dict()
    handle = open(idmappingfile,'rb')
    handle.readline()
    for line in handle:
        fields = line.strip().split('\t')
        uniprotdict[fields[0]]=fields[1]
    handle.close()
    print('Uniprot Id Mapping read for %s.\n' % taxon )
    return(uniprotdict)
    

def __write_uniprot_id__(taxon):
    #This function outputs a file with the list of accessions 
    #to put through the web tool   
    #Only need to run one time!
    if taxon in ['208963','7227','237561']:
        filename = 'mapping.'+str(taxon)+'.map'
    else:
        filename = 'sp_species.'+str(taxon)+'.map'
    handle = open('./Mapping files/'+filename,'r')
    if taxon in file_mapping_available:
        pass
    else:
        outfile = open('./uniprotID_'+str(taxon)+'.txt','w')
        for line in handle:
            fields = line.strip().split('\t')
            outfile.write('%s\n' % fields[1])
    outfile.close()
    

    
def read_target_mapping(taxon):
    targetdict = dict()
    if taxon in ['208963','7227','237561']:
        filename = 'mapping.'+str(taxon)+'.map'
    else:
        filename = 'sp_species.'+str(taxon)+'.map'
    handle = open('./Mapping files/'+filename,'r')

    for line in handle:
        fields = line.strip().split('\t')
        targetdict[fields[1]]=fields[0]
    handle.close()
    return(targetdict)


def type_converter(oldtype):
    if oldtype=='NK':
        newtype='type1'
    elif oldtype=='LK':
        newtype = 'type2'
    elif oldtype=='ALL':
        newtype = 'all'
    return(newtype)
                
def mapping(benchmarkFolder,preBench,uniprotdict,taxon,Type,onto):
    #benchmarkFolder is the folder where newCAFA assessment software can take
    #preBench is a benchmark file generated by ataur's code (or Huy's code)
    #id_mapping is a uniprot mapping file
    #target_mapping is a CAFA3 mapping file
    #This function maps the Ataur (or Huy) Benchmark results into newCAFA assessment software formats
    targetdict = read_target_mapping(taxon)
    groundtruth = open(benchmarkFolder+'/groundtruth/leafonly_'+Type+'_'+taxon+'_'+onto.upper()+'.txt','w')
    lists = open(benchmarkFolder+'/lists/'+onto+'_'+taxon_name_converter(taxon)+'_'+type_converter(Type)+'.txt','w')
    if os.path.getsize(preBench)==0:
        print('No Benchmark Found in %s_%s_%s\n' % (taxon,Type,onto))
    else:
        ahandle = open(preBench,'r')
        cafaIDs = set()
        for line in ahandle:
            fields = line.strip().split('\t')
            try:
                uniprotID = uniprotdict[fields[0]]
            except KeyError:
                print('Benchmark not in uniprot mapping: %s\n' % fields[0])
                continue
            try:
                cafaID = targetdict[uniprotID]
                groundtruth.write('%s\t%s\n' % (cafaID,fields[1]))
                cafaIDs.add(cafaID)
            except KeyError:
                print('Benchmark not in target file:%s.\n' % uniprotID)
        for cid in cafaIDs:
            lists.write('%s\n' % cid)
        ahandle.close()
    groundtruth.close()
    lists.close()
        
    
def run(benchmarkFolder,taxon_id,f,Type,onto):
    if taxon_id not in file_mapping_available:
        print('Mapping file of %s not directly available from uniprot.\n' % taxon_id)
        uniprotdict = read_uniprot_IdMapping_fromWeb('./uniprot_ac_to_id_%s.tab' % taxon_id,taxon_id)
        mapping(benchmarkFolder,f,uniprotdict,taxon_id,Type,onto) 
    else:
        uniprotdict=read_uniprot_IdMapping('./uniprot_ac_to_id_%s.map' % taxon_id,taxon_id)
        mapping(benchmarkFolder,f,uniprotdict,taxon_id,Type,onto)         
        
        
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Map UniProt Accession ID to CAFA ID', )
    parser.add_argument('file',help='Input prebenchmark file')
    parser.add_argument('-b',dest='benchFolder',help='Final Benchmark Folder',required=True)
    parser.add_argument('-s',dest='species',help='Input list of species, in NCBI taxonomy ID',required=True)
    parser.add_argument('-t','--t',dest='type',help = 'Input evaluation type: No Knowledge or Limited Knowledge', choices=['NK','LK'],required=True)    
    parser.add_argument('-o','--o',dest='onto',help='Input ontology',choices=['bpo','cco','mfo'],required=True)
    args = parser.parse_args()
    run(args.benchFolder,args.species,args.file,args.type,args.onto)
    