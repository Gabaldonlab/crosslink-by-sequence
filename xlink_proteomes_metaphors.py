#!/usr/bin/env python
desc="""
This script is based on original Leszeck script to cross link proteomes src/xlink_proteomes.py
The script takes as input a reference proteome file with protIds in the fasta headers and proteome for cross linking. 
The output of the script is a ready to load table with ext2meta linked protomes, fasta file with re-named sequences, and fasta file with orphan sequences.
"""
epilog="""Author:
l.p.pryszcz@gmail.com
Barcelona/Borki, 14/11/2012
Barcelona, 5/09/2013

major changes by A.Vlasova
Barcelona, 17/02/2021 
"""
#TODO - add parameter to choose to run blat or diamond. 
#  
import argparse, os, sys
import gzip, hashlib, time
import subprocess as sp

from Bio        import SeqIO
from datetime   import datetime
#from threading  import activeCount, Thread
#from tempfile import NamedTemporaryFile
from multiprocessing import Pool
import pandas as pd
import numpy as np 

def get_md5_fasta(fname):
    """md52prot 
    """
    hash2prot  = {} # {  FASTA_md5: [ protid1, protid2 ]  ...  }
    with gzip.open(fname, mode='rt') as handle:
        for r in SeqIO.parse(handle, 'fasta'):
            #add prot2fasta
            extID = r.id
            seq   = str(r.seq)
            #add hash2prot
            md5   = hashlib.md5(seq.encode('utf-8')).hexdigest()
            if md5 not in hash2prot:
                hash2prot[md5] = []
            hash2prot[md5].append(extID)       
    print(dict(list(hash2prot.items())[0:2]))
    return hash2prot

def makedb_diamond(tmpDir, refFile):           
    outputFileName = tmpDir + os.path.basename(refFile)+'db'    
    cmd = "diamond makedb -d %s --in %s " % ( outputFileName,refFile )    
    #print(cmd)
    if not os.path.isfile(outputFileName+'.dmnd'):
        run_command(cmd, False)
        print(cmd)
    else:
        print("[INFO] Diamond DB file exists for ", refFile)

def get_sequence_length(fileName):
    seqLength  = {} # {  FASTA_md5: [ protid1, protid2 ]  ...  }
    with gzip.open(fileName, mode='rt') as handle:
        for r in SeqIO.parse(handle, 'fasta'):
            #add prot2fasta
            extID = r.id
            seqL   = len(str(r.seq))
            seqLength[extID]=seqL
    return seqLength

def xLink_diamond(ext2meta, tmpDir, refFile, missingFn, covTh, idTh, threads):
    
    #run diamond
    dbFile = tmpDir + os.path.basename(refFile)+'db'       
    outputFileName = tmpDir + os.path.basename(missingFn) + 'diamond'   
    
    #-tileSize=8 caused `Internal error genoFind.c 2225` error
    cmd = "diamond blastp --db %s --query %s --out %s --outfmt 6 --header --strand both --threads %s" % ( dbFile,missingFn,outputFileName, threads ) 
    #print(cmd)
    #TODO - check existance of the file
    run_command(cmd, False)
    #get length of each sequence in missing file 
    lenSequence = get_sequence_length(missingFn)
    lenSequenceSubj = get_sequence_length(refFile)

    
    #parse diamond output
    q2t={}
    Qset, Tset = set(), set()
    for line in open(outputFileName):
        line = line.strip()
        if not line or line.startswith("#"):
            continue        
        #parse with error capturing
        #Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
        lData = line.split('\t')
        if len(lData) == 12:
            Q,T,identity,aLength,nMis,nGaps,sQ,eQ,sS,eS,expV,score = lData
            aLength = int(aLength)
            identity=float(identity)/100
            if lenSequence[Q] > lenSequenceSubj[T]:
                length = lenSequence[Q]
            else:
                length = lenSequenceSubj[T]
            cov = 1.0*(aLength) / length
            #store Q and T
            Qset.add(Q)
            Tset.add(T)
        #skip lines with errors
        else:
            if verbose: 
                sys.stderr.write("  Error: blat2xlinks: Expected 21 elements in line, got %s %s!\n" % (len(lData), str(lData)))
            continue 
    
        #filter low identity or low coverage hits
        if cov < covTh or identity < idTh: 
            continue
    
        Qprotid = Q 
        Tprotid = T
        iScore = (identity+cov) / 2
        Tdata  = (Tprotid, iScore)       
       
        #save only the best match for each query, more than one allowed, if best matches has got the same score
        if Qprotid in q2t:
            #append matches if same iScore
            if   iScore == q2t[Qprotid][0][1]: 
                #add entry in t2q
              #  if Tprotid not in t2q:
              #      t2q[Tprotid]=[]
              #  t2q[Tprotid].append((Qprotid, iScore))
                #add entry to q2t
                q2t[Qprotid].append(Tdata)
            #update matches if better iScore
            elif iScore  > q2t[Qprotid][0][1]: 
                q2t[Qprotid] = [Tdata]                                
        else:
            q2t[Qprotid] = [Tdata]
            #add entry in t2q
            #if Tprotid not in t2q:
            #    t2q[Tprotid]=[]
            #t2q[Tprotid].append((Qprotid, iScore))
    #print("T2Q:", len(t2q))
    #print(dict(list(t2q.items())[0:2]))
    #print("Q2T:")
    #print(dict(list(q2t.items())[0:2]))
    
    df = pd.DataFrame.from_dict(q2t, orient ="index")
    #print(df)
    if len(df.index)>1:
        df=df.unstack().dropna().sort_index(level=1)
        df=df.reset_index(level=[0])
        df["extId"]=df.index
        df=df.reset_index()

        df.columns =["index","oldIndex","value","extId"]

        ext2metaBlat=pd.DataFrame(df['value'].tolist())
        ext2metaBlat.columns=["protid","score"]
        ext2metaBlat["extid"]=df["extId"]
        #extid score  protid
        ext2metaBlat=ext2metaBlat[["extid", "score",  "protid"]]
        #print(ext2metaBlat.head())
        #print(ext2meta.head())
        ##select only rows that were not matched by md5 
        ext2metaBlat=ext2metaBlat[~ext2metaBlat["extid"].isin(ext2meta["extid"])]
        ext2meta=pd.concat([ext2meta, ext2metaBlat], axis=0)
    return ext2meta


def xLink_blat(ext2meta, tmpDir, refFile, missingFn, covTh, idTh):
    
    #run blat   
    outputFileName = tmpDir + os.path.basename(missingFn)    
    psl = "%s.psl" % outputFileName
    #-tileSize=8 caused `Internal error genoFind.c 2225` error
    cmd = "blat -prot -out=psl -minIdentity=90 -minMatch=1 -out=psl -noHead %s %s %s" % ( refFile,missingFn,psl )    
    # print(cmd)
    #TODO - check existance of the file
    run_command(cmd, False)
    #parse blat output
    q2t={}
    Qset, Tset = set(), set()
    for line in open(psl):
        line = line.strip()
        if not line:
            continue        
        #parse with error capturing
        lData = line.split('\t')
        if len(lData) == 21:
            match,mMatch,rMatch,Ns,Qgapcount,Qgapbases,Tgapcount,Tgapbases,strand,Q,Qsize,Qstart,Qend,T,Tsize,Tstart,Tend,bcount,bsizes,qStarts,tStarts = lData
            match,mMatch,Qsize,Tsize = int(match),int(mMatch),int(Qsize),int(Tsize)
            # take longer query or target size as general size of compared proteins
            if Qsize>Tsize:
                length = Qsize
            else:
                length = Tsize
            #length=Q_size
            cov = 1.0*(match+mMatch) / length
            identity = 1.0*match / length
            #store Q and T
            Qset.add(Q)
            Tset.add(T)
        #skip lines with errors
        else:
            if verbose: 
                sys.stderr.write("  Error: blat2xlinks: Expected 21 elements in line, got %s %s!\n" % (len(lData), str(lData)))
            continue 
    
        #filter low identity or low coverage hits
        if cov < covTh or identity < idTh: 
            continue
    
        Qprotid = Q 
        Tprotid = T
        iScore = (identity+cov) / 2
        Tdata  = (Tprotid, iScore)   
       
        #save only the best match for each query, more than one allowed, if best matches has got the same score
        if Qprotid in q2t:
            #append matches if same iScore
            if   iScore == q2t[Qprotid][0][1]: 
                #add entry in t2q
              #  if Tprotid not in t2q:
              #      t2q[Tprotid]=[]
              #  t2q[Tprotid].append((Qprotid, iScore))
                #add entry to q2t
                q2t[Qprotid].append(Tdata)
            #update matches if better iScore
            elif iScore  > q2t[Qprotid][0][1]: 
                q2t[Qprotid] = [Tdata]                                
        else:
            q2t[Qprotid] = [Tdata]
            #add entry in t2q
            #if Tprotid not in t2q:
            #    t2q[Tprotid]=[]
            #t2q[Tprotid].append((Qprotid, iScore))
    #print("T2Q:", len(t2q))
    #print(dict(list(t2q.items())[0:2]))
    #print("Q2T:")
    #print(dict(list(q2t.items())[0:2]))
    
    df = pd.DataFrame.from_dict(q2t, orient ="index")
    df=df.unstack().dropna().sort_index(level=1)
    df=df.reset_index(level=[0])
    df["extId"]=df.index
    df=df.reset_index()

    df.columns =["index","oldIndex","value","extId"]

    ext2metaBlat=pd.DataFrame(df['value'].tolist())
    ext2metaBlat.columns=["protid","score"]
    ext2metaBlat["extid"]=df["extId"]
    #select only rows that were not matched by md5 
    ext2metaBlat=ext2metaBlat[~ext2metaBlat["extid"].isin(ext2meta["extid"])]
    ext2meta=pd.concat([ext2meta, ext2metaBlat], axis=0)
    #ext2meta2 = pd.DataFrame.from_dict(q2t)  
    return ext2meta

def xLink_md5(hash2meta, hash2ext):
    """Link proteins using md5. 
    Return ext2meta for given dbID
    """
    ext2meta =  pd.DataFrame(columns=['extid','score','protid'])  
    for md5 in hash2ext:
        #if same seq already present
        if md5 in hash2meta: 
            #metaID is a list, can be a few elements
            metaID = hash2meta[md5]            
            for extID in hash2ext[md5]:
                #add ext2meta entry    
                for metIDEl in metaID:                    
                    ext2meta=ext2meta.append(pd.DataFrame([[extID,'1',metIDEl]], columns=ext2meta.columns))
                   
    #print(dict(list(ext2meta.items())[0:2]))
    #print(ext2meta.head())
    return ext2meta


#outDir, tmpDir, refFile, fnames, taxid, minCoverage, minIdentity, protLimit
def process_taxid(outDir, tmpDir, refFile, fileName,hash2meta, taxid, minCoverage, minIdentity, protLimit):    
    
    name = os.path.basename(fileName)
    print("[INFO] Xlink ", name)
    dbid=0
    version=0
    # dbid, strtaxa, prefix, _ = name.split(".", 3)
    # if prefix == "longest" or prefix == "faa":
    #     prefix = 0
    # version = 0 
    # if dbid == '8':
    #     version = prefix
    #     prefix = prefix +"."+ strtaxa

    # for a few DBs we may have a few strains for the same species, that does not have versions 
    # (in March 2021 I've observed it for eggnog, treefam, and evoclust)
    #this lead that result files are saved under same name and the last process rewrite output of the previous processes 
    # Example:
    #/gpfs/projects/bsc40/project/pipelines/metaphors/metaphors-db-2019/download/eggnog4.5/fasta.v5/10.521097.faa.gz, 
    #/gpfs/projects/bsc40/project/pipelines/metaphors/metaphors-db-2019/download/eggnog4.5/fasta.v5/10.873517.faa.gz,
    #and refernce file /gpfs/projects/bsc40//project/pipelines/metaphors/metaphors-db-2019/data/fasta.reference//10.1018.faa.gz
    #belong to the same taxa 1018

   # if dbid in [ '10', '12', '13' ]:
   #     prefix = strtaxa 
    
    #filePrefix=str(dbid)+"."+str(taxid)+"."+str(prefix)
    filePrefix=name
      
    hash2ext =  get_md5_fasta(fileName)
    
    ext2metaTmp = xLink_md5(hash2meta, hash2ext)

    numberAll = run_command_with_return("zcat %s |grep -c '>'"%fileName)
    numberAll = numberAll[0].decode('utf-8').strip()
    numberAll = int(numberAll)
    numberMatched = ext2metaTmp["extid"].nunique()
    
    rest = numberAll - numberMatched
    print("Rest ", rest)

    if rest>0:
        #run blat comparison between two files in the temporary folder 
        #this function returns ext2meta dictionary with keys - protein names from fileName and values - protein names from refFile 
        #ext2metaTmp = xLink_blat(ext2metaTmp, tmpDir, refFile, fileName, minCoverage, minIdentity)        
        ext2metaTmp = xLink_diamond(ext2metaTmp, tmpDir, refFile, fileName, minCoverage, minIdentity,1)        
    

    ext2metaTmp["dbid"]=dbid
    ext2metaTmp["version"]=version
    
    ext2metaTmp=ext2metaTmp[["dbid","extid","version","protid","score"]]

    #format ext2meta to the format compartible with the current ext2meta schema in the DB
    ext2metaFn  = os.path.join(outDir, '%s.ext2meta.tbl.gz'%filePrefix)
    ext2metaTmp.to_csv(ext2metaFn, compression="gzip", index=False, sep="\t", header=False)
    #format fasta file with orphans
    #outOrphans = outDir+filePrefix+".orphans.faa.gz"
    #print("[INFO] writing orphans for ",name )
    #saveGzipFasta(fileName,ext2metaTmp["extid"].unique(), outOrphans, 'invert' )    

    #at the same time calculate number of matched and non-matched proteins    
    #print some stats with % of orphans per file 
    numberMatched = ext2metaTmp["extid"].nunique()
    numberOrphans = numberAll - numberMatched 
    percOrphans = numberOrphans *100 /numberAll

    logFn = os.path.join(outDir,'%s.log' % filePrefix)
    log   = open(logFn,'wt')
    print("%s\t%s\t%s\t%s\t%.2f"%(filePrefix, numberAll, numberMatched,numberOrphans, percOrphans), file=log)    
    log.close()
    
    return filePrefix

def get_taxids(fileName):
    """Return taxids and species2strains"""
    taxid = fileName.split('/')[-1].split('.')[1]
    return taxid

def fasta_generator(outDir, tmpDir, refFile, fnames, taxid,  minCoverage, minIdentity, protLimit):
    for fastaFile in fnames:
        yield(outDir, tmpDir, refFile, fastaFile, taxid,  minCoverage, minIdentity, protLimit)

def saveGzipFasta(fileName, idList, outName, param):
    handle_in = gzip.open(fileName, "rt")
    handle_out = gzip.open(outName, "wt")    
    for record in SeqIO.parse(handle_in, "fasta"):            
        if np.isin(record.id, idList):
            if param == 'invert':
                continue
            else:               
                handle_out.write(record.format("fasta"))
        else:
            if param == 'invert':
                handle_out.write(record.format("fasta"))
            else:
                continue

            
        
    
def process(fnames, refFile, outDir, tmpDir, nprocs, minCoverage, minIdentity, protLimit):
    """xlink proteomes """ 
    ###get taxids and species2strains
    #taxid = get_taxids(refFile)    
    taxid=99999999
    hash2meta = get_md5_fasta(refFile)
    ###process files 
    #print("[INFO] Processing taxa ", taxid)
    makedb_diamond(tmpDir, refFile)

    p = Pool(nprocs)   
    #tdata = fasta_generator(outDir, tmpDir, refFile, fnames, taxid, minCoverage, minIdentity, protLimit)
    tdata=[]
    for fastaFile in fnames:
        tdata.append([outDir, tmpDir, refFile, fastaFile,hash2meta, taxid,  minCoverage, minIdentity, protLimit])
                           
    for i, data in enumerate(p.starmap(process_taxid, tdata), 1):
        if not data:
            continue
        fileNamePrefix = data
        sys.stderr.write(" %s / %s  %s    \r"%(i, len(fnames), fileNamePrefix))
        print("[INFO] done %s "%(fileNamePrefix))     

def create_folder(name):    
    if not os.path.exists(name):      
        cmd = "mkdir "+name
        try:
            run_command(cmd,False)
        except:
            print("Unable to create directory ",name)

def run_command(cmd,ommit):
    if ommit:
        try: process = sp.Popen(cmd,shell=True)
        except: pass
        process.communicate("Y\n")
        if process.wait() != 0: print("Error ocurred, but you chose to ommit it")
    else:
        try: process = sp.Popen(cmd,shell=True)
        except OSError as e: sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0: sys.exit("ERROR: Execution cmd failed")
        
def run_command_with_return(cmd):
    try: process = sp.Popen(cmd,shell=True,stdout = sp.PIPE).stdout
    except: sys.exit()
    return process.readlines()

def main():

    usage  = "%(prog)s [options] -f */fasta/*.faa.gz"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", dest="verbose",     default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.2')
    #parser.add_argument("-f", dest="fastas",      nargs="+", #type=file,
    #                    help="proteome fasta files")
    parser.add_argument("-f", dest="fastas",   nargs="+", required=True,  #type=file,
                        help="proteome fasta files")                    
    #parser.add_argument("-f", "--fastaFile", dest="fastaFile", required=True, default="",
    #                    help="Tabular separated file with spTaxid, Taxid, path, and some other info")
    parser.add_argument("--refSpecie", dest="refSpecie",   required=True, 
                        help="Reference species fasta file")
    parser.add_argument("-o", dest="outDir",      default='fasta.crosslinked', type=str,
                        help="output directory             [%(default)s]")  
    parser.add_argument("--tmp", dest="tmpDir",      default='', type=str,
                        help="TMP directory to store blat psl files [%(default)s]")                            
    parser.add_argument("-c", dest="minCoverage", default=0.95, type=float,
                        help="min. coverage for blat hits  [%(default)s]")
    parser.add_argument("-i", dest="minIdentity", default=0.98, type=float,
                        help="min. identity for blat hits  [%(default)s]")
    parser.add_argument("-l", dest="protLimit",   default=0, type=int,
                        help="only taxa having >l proteins [%(default)s]")
    parser.add_argument("-t", "--threads",        default=4, type=int,
                        help="number of cores [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
  
    #create out and tmp directory if they dont exists
    create_folder(o.outDir)
    create_folder(o.tmpDir)
  
    #xLink proteomes   
    process(o.fastas,o.refSpecie, o.outDir, o.tmpDir, o.threads, o.minCoverage, o.minIdentity, o.protLimit) 
            
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )

