#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Pedro Navarro$
$Authors: Pedro Navarro$
--------------------------------------------------------------------------
"""
from __future__ import print_function

import sys
import os
import csv
import getopt 
import glob
import re
from configobj     import ConfigObj
#DEPRECATED: from cluster     import HierarchicalClustering
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.neighbors import kneighbors_graph
    
import     msproteomicstoolslib.format.speclib_db_lib         as         speclib_db_lib

def usage() :
    print("")
    print("spectrast_cluster.py")
    print ("-" * 20)
    print("This script clusters by iRT peptide spectra from spectraST, and writes the clusters in separate .sptxt files")
    print("")
    print("Usage: ")
    print("python cluster_spectrast.py [options] spectrast_file(s)")
    print("-h            --help        Display this help")
    print("-d    distance        --distance        Cluster distance. Default : 1")
    print("")

def clusterRT_ward(values) :
    if len(values) == 0 : return []

    v = sorted([[val] for val in values])
    
    #connectivity = kneighbors_graph(np.asarray(v), n_neighbors=3)
    ward = AgglomerativeClustering(n_clusters=2).fit(np.asarray(v))
    labels = ward.labels_

    curr_l = -2
    cl_output = []
    curr_cluster = []
    for i,l in enumerate(labels) :
        if l != curr_l :
            if len(curr_cluster) > 0 : cl_output.append(curr_cluster)
            curr_l = l
            curr_cluster = []
        curr_cluster.append(values[i])
    cl_output.append(curr_cluster)
    
    return cl_output

def clusterRT_meanshift(values,bandwidth) :
    if len(values) == 0 : return []

    v = sorted([[val] for val in values])
    
    bandwidth = estimate_bandwidth(np.asarray(v), quantile=0.1, n_samples=int(len(v)))
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, min_bin_freq = 2, cluster_all = False)
    ms.fit(np.asarray(v))
    labels = ms.labels_

    curr_l = -2
    cl_output = []
    curr_cluster = []
    for i,l in enumerate(labels) :
        if l != curr_l :
            if len(curr_cluster) > 0 : cl_output.append(curr_cluster)
            curr_l = l
            curr_cluster = []
        curr_cluster.append(values[i])
    cl_output.append(curr_cluster)
    
    return cl_output
    

def clusterRT(values, rt_maximal_distance, algorithm=False) :
    if algorithm == "ward" : return clusterRT_ward(values)
    if algorithm == "meanshift" : return clusterRT_meanshift(values, rt_maximal_distance)
     
    #If not any sample is over the threshold --> return the list as it is
    if len(values) == 0 : return []
    
    ## DEPRECATED : we don't use anymore the cluster library
    #cl = HierarchicalClustering(values, lambda x,y: abs(x-y))
    #cl_output = cl.getlevel(rt_maximal_distance)     # get clusters of items closer than rt_maximal_distance

    values = sorted(values)
    v = [[val] for val in values]
    #print v
    db = DBSCAN(eps=rt_maximal_distance,min_samples=2).fit(np.asarray(v))
    labels = db.labels_
    #print labels
    
    curr_l = -2
    cl_output = []
    curr_cluster = []
    for i,l in enumerate(labels) :
        if l != curr_l :
            if len(curr_cluster) > 0 : cl_output.append(curr_cluster)
            curr_l = l
            curr_cluster = []
        curr_cluster.append(values[i])
    cl_output.append(curr_cluster)
    
    #print cl_output
    return cl_output




def main(argv) :

    distance        = 1.0
    algorithm       = False
    
    #Get options
    try:
        opts, args = getopt.getopt(argv, "hd:i:t:a:",["help","distance","algorithm"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    argsUsed = 0
    for opt,arg in opts:
        if opt in ("-h","--help") :
            usage()
            sys.exit()
        if opt in ("-d","--distance") :
            distance = float(arg)
            argsUsed += 2
        if opt in ("-a","--algorithm") :
            algorithm = arg
            argsUsed += 2
    
    
    sptxtfiles_pat = argv[argsUsed:]
    sptxtfiles = []
    for pat in sptxtfiles_pat :
        sptxtf = glob.glob(pat)
        for file in sptxtf : sptxtfiles.append(file)


    for sptxtfile in sptxtfiles :
        transitions = []
        print("Reading : " , sptxtfile)
        assert sptxtfile[-6:] == '.sptxt'
        if not os.path.exists(sptxtfile):
            print("The file: %s does not exist!" % sptxtfile)
            sys.exit(2)
            
        library_key = 99
        spectrastlib = speclib_db_lib.Library(library_key)

        num_spectrum = 0
        offset = spectrastlib.get_first_offset(sptxtfile)
        last_offset = -100

        #Get all the peptide sequences and retention times to cluster them. Keep the spectrum number associated. 
        
        peptide_spectra = {}  # { "SEQUEN[Pho]CE" : {last_offset1 : RT1  , last_offset2 : RT2 } , "SEQUENCE" : {last_offset3 : RT3 , last_offset4 : RT4 , ... } , ...  }
        
        
        while ( offset - last_offset > 10) :
            last_offset = offset
            offset , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)

            #for property, value in vars(spectrum).iteritems():
            #    if property    in ['compress_spectra' ] : continue
            #    print property, ": ", value
            #sys.exit()

            sequence     = spectrum.name.split('/')[0]
            z_parent     = float(spectrum.name.split('/')[1])
            if spectrum.RetTime_detected:
              rt             = spectrum.RetTime
            if spectrum.iRT_detected:
              rt             = spectrum.iRT
            if not spectrum.RetTime_detected and not spectrum.iRT_detected:
              print("No RT/iRT was detected for %s" % spectrum.name)
              sys.exit(2)
            
            if sequence in list(peptide_spectra.keys()) :
                peptide_spectra[sequence][last_offset] = rt
            else :
                peptide_spectra[sequence] = { last_offset : rt }
            
        
        max_num_of_clusters = 0
        peptide_spectra_cl = {}
        
        print("cluster spectra by iRTs...")
        for sequence, spectra in peptide_spectra.items() :
            print(sequence)  #, spectra
            rt_clusters = clusterRT(list(spectra.values()), distance, algorithm = algorithm)
            
            if len(rt_clusters) > max_num_of_clusters : max_num_of_clusters = len(rt_clusters)
            
            peptide_spectra_cl[sequence] = {}
            
            for spectrum, rt in spectra.items() :
                # Determine cluster number for this rt
                cl_index = -1
                for index, cluster in enumerate(rt_clusters) :
                    #print index, rt, cluster
                    cl = cluster
                    if not isinstance(cluster,list) : cl = [cluster] 
                    if rt in cl : cl_index = index
                
                #store cluster index in a dictionary      
                peptide_spectra_cl[sequence][spectrum] = cl_index
        
            
        splitfiles = [ open(sptxtfile[:-6]+"_"+str(x+1)+".sptxt",'w') for x in  range(max_num_of_clusters)  ]        

        #init the files by using the original header
        print("%s files will be created." % max_num_of_clusters)
        
        original_header = spectrastlib.get_fileheader(sptxtfile)
        for file in splitfiles :
            for line in original_header : file.write(line)
                     
        for sequence, spectra in peptide_spectra_cl.items() :
            for spectrum in spectra :
                sp = spectrastlib.get_rawspectrum_with_offset(sptxtfile,spectrum)  #get the spectrum
                
                for line in sp :
                    #Add suffix to the protein name
                    if 'Comment:' in line[:8] : 
                        line_bcp = line
                        mm = re.search( 'Protein=(.*?)\s', line )
                        if not mm: break 
                        split_idx = line.index('Protein=') + 8
                        line_before_split = line[:split_idx]
                        line_after_split = line[split_idx:]
                        line = line_before_split + "Subgroup_" + str(peptide_spectra_cl[sequence][spectrum]) + "_" + line_after_split
                        
                    splitfiles[peptide_spectra_cl[sequence][spectrum]].write(line)
                    
        for file in splitfiles :
            file.close()
        
        print("done.")
        
if __name__ == '__main__':
    main(sys.argv[1:])
