#!/usr/bin/env python
"""Fetch all entries for vertebrate ND3s."""

import sys
from Bio      import Entrez
from datetime import datetime
import os

# please change to a valid email 
# Entrez.email = "xxx@xxx.com"
query="\"NADH dehydrogenase subunit 3\" AND Vertebrata[Organism] AND mitochondrion[filter]"
db="nuccore"
retmax=10**9
retmode='text'
rettype='gb'
batchSize=100

#get list of entries for given query
sys.stderr.write( "Getting list of GIs for term=%s ...\n" % query )
handle = Entrez.esearch( db=db,term=query,retmax=retmax )
giList = Entrez.read(handle)['IdList']

#print info about number of proteins
sys.stderr.write( "Downloading %s entries from NCBI %s database in batches of %s entries...\n" % ( len(giList),db,batchSize ) )

#post NCBI query
search_handle     = Entrez.epost( db, id=",".join( giList ) )
search_results    = Entrez.read( search_handle )
webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 

#fecth all results in batch of batchSize entries at once
for start in range( 0,len(giList),batchSize ):
    #print info
    tnow = datetime.now()
    sys.stderr.write( "\t%s\t%s / %s\n" % ( datetime.ctime(tnow),start,len(giList) ) )
    out_f = "out_"+str(start)+'.gb'
    if not os.path.isfile(out_f):
        try:
            handle = Entrez.efetch( db=db,retmode=retmode,rettype=rettype,retstart=start,retmax=batchSize,webenv=webenv,query_key=query_key )
            with open(out_f, "w") as out:
                out.write( handle.read() )
        except:
            sys.stderr.write("error download {} + {} sequences".format(start,batchSize))
            

