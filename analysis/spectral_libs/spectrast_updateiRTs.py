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
from __future__ import division
from __future__ import print_function

import sys, csv
import os
import getopt, glob
from configobj import ConfigObj

from msproteomicstoolslib.math.LinearRegression import SimpleLinearRegression

def usage() :
    print("")
    print("spectrast_updateiRTs.py")
    print ("-" * 23)
    print("This script updates the RetTime values of a spectraST files with different values (i.e. iRTs) calculated with a linear model.")
    print("")
    print("Usage: ")
    print("python spectrast_updateiRTs.py [options] spectrast_file(s)")
    print("-a            --auto-align  Calculates automatically the alignment models based on the internal identifications of the calibration peptides")
    print("-h            --help        Display this help")
    print("-i    iRT_models_file    --irtmodel    File with the iRT models")
    print("-p    iRT_peptides_file    --irtpeptides File containing the peptides and their iRTs to do the linear model alignment.")
    print("-t   time-scale            Options: minutes, seconds. Default: seconds.")
    print("")

def read_irtmodels(file, useMinutes = False) :
    irtmodels = {}  # { file1 : (a1,b1) , file2 : (a2,b2), ...  }
    
    fs = csv.reader(open(file), delimiter="\t")

    for cnt, srow in enumerate(fs):
        if '#' in srow[0]:
            continue
        if not useMinutes:
            irtmodels[srow[0]] = (float(srow[1]), float(srow[2]))
        else:
            irtmodels[srow[0]] = (float(srow[1]), srow[2] / 60.0)
    
    print("iRT models : ", irtmodels)
    
    return irtmodels

def cal_irtmodels(irts_in_samples) :
    '''
    This calculates a linear regression retention time model for each sample. 
    
    The input should be :
    irts_in_samples = { 'Sample1' : { 'Pep1' : (iRT, RT,Int) , 'Pep2' : (iRT, RT,Int), ... } , ...  }
    
    The output is :
    irtmodels = {}  # { file1 : (b1,a1) , file2 : (b2,a2), ...  }  , where a == intersection, b == slope
    '''   
    
    R2threshold     = 0.98
    max_outliers    = 2
    
    irtmodels       = {}  # { file1 : (a1,b1) , file2 : (a2,b2), ...  }
    
    for sample in irts_in_samples : 
        data = [(RT, iRT) for (iRT,RT,Int) in irts_in_samples[sample].values()]  # Retrieve the data in the format [ (x1,y1) , (x2,y2) , ... ]

        linRegr = SimpleLinearRegression(data)
        if not linRegr.run():
            print("...error: failed to calculate parameters")
            return

        print("...linear function of %s is then %s . (r** = %f)" % (  sample,linRegr , linRegr.r**2) )
        
        if linRegr.r**2 < R2threshold :
            #Remove the farthest outliers (a maximum of max_outliers)
            curr_outlier = (0.0,0.0)
            for ol in range(0,max_outliers) :
                max_dist = 0.0
                for (x, y) in data :
                    yp = linRegr.function(x)
                    dist = abs(yp - y)
                    if dist >= max_dist :
                        max_dist = dist
                        curr_outlier = (x, y)
                print("Removing outlier (%f, %f)..." % ( x , y ))
                data.remove(curr_outlier)
                linRegr = SimpleLinearRegression(data)
                if not linRegr.run():
                    print("...error: failed to calculate parameters")
                    return
                print("... new linear function of %s is then %s . (r** = %f)" % (  sample,linRegr , linRegr.r**2) )
                if linRegr.r**2 >= R2threshold : break
                    
        irtmodels[sample] = ( linRegr.a , linRegr.b )
     
    return irtmodels

def search_irtPeptides(file, irtPeptides) :
    import re
    
    irts_in_samples = {}  # This will contain : { 'Sample1' : { 'Pep1' : (iRT, RT,Int) , 'Pep2' : (iRT, RT,Int), ... } , ...  }
    
    f = open(file,'r')
    
    curr_peptide = ''
    curr_sample  = ''
    curr_RetTime = 0.0
    curr_Intensity = 0.0
    
    for row in f :
        if 'Name:' in row[:5] :
            curr_peptide = row[6:].strip()
        if curr_peptide not in irtPeptides : continue
        
        if 'Comment:' in row[:8] :
            mm = re.search( 'RawSpectrum=(.*?)\.', row )
            if not mm : continue
            curr_sample = mm.group(1)
            
            mm = re.search( 'RetentionTime=(.*?)\s', row )
            
            if not mm : continue
            
            if ',' in mm.group(1)[:] :
                mm_split        = mm.group(1).split(',')
                mm_split_float  = [float(f) for f in mm_split]
                curr_RetTime    = mm_split_float[1] # I pick the RT in the middle (average?)
            else : curr_RetTime  = float( mm.group(1) )
            
            mm = re.search( 'OrigMaxIntensity=(.*?)\s', row )
            
            if mm : curr_Intensity = float(mm.group(1))
            else  : curr_Intensity = 0.0
            
            if curr_sample not in irts_in_samples : irts_in_samples[curr_sample] = { curr_peptide : ( irtPeptides[curr_peptide] , curr_RetTime , curr_Intensity) }
            elif curr_peptide not in irts_in_samples[curr_sample] : irts_in_samples[curr_sample][curr_peptide] = ( irtPeptides[curr_peptide] , curr_RetTime , curr_Intensity )
            else : 
                (irt, rt, prev_intensity) = irts_in_samples[curr_sample][curr_peptide]
                if curr_Intensity >= prev_intensity : irts_in_samples[curr_sample][curr_peptide] = ( irtPeptides[curr_peptide] , curr_RetTime , curr_Intensity )
    
    return irts_in_samples
            
def update_iRTs(file, models) :
    import re
    
    f = open(file,'r')
    f_new = open(file + '.tmp', 'w')

    for row in f :
        if 'Comment:' not in row[:8] : 
            f_new.write(row)
            continue
        
        #Get which model to use
        modelkey = ''
        curr_model =  ( 0.0 , 0.0 )
        for key, model in models.items() :
            #print key, model
            keym = 'RawSpectrum=' + key + '.'
            if keym in row : curr_model = model
        
        RetTime = 0.0 
        mm = re.search( 'RetentionTime=(.*?)\s', row )
        
        if mm:
            if ',' in mm.group(1)[:] :
                mm_split        = mm.group(1).split(',')
                mm_split_float  = [float(f) for f in mm_split]
                RetTime         = mm_split_float
            else : RetTime  = float( mm.group(1) )
            
            if isinstance(RetTime,list) : 
                for idx, rt in enumerate(RetTime) : RetTime[idx] = curr_model[0] + RetTime[idx] * curr_model[1]
            else : RetTime = curr_model[0] + RetTime * curr_model[1]
        
            #Split the row into what it is BEFORE the RetTime, and what it is AFTER the RetTime
            split_idx = row.index('RetentionTime=')
            row_before = row[:split_idx]
            split_idx = split_idx + len(mm.group(0))
            row_after = row[split_idx:]
            
            row = row_before + 'RetentionTime='
            if isinstance(RetTime, list) :
                row = row + ','.join([str(item) for item in RetTime]) 
            else : row = row + RetTime
            row = row + ' ' + row_after
        
        f_new.write(row)
        
    f_new.close()
    
    base = os.path.splitext(file)[0]
    os.rename(file + '.tmp', base + "_iRT.sptxt")
    
    return


def read_irtPeptides(filename) :
    """Reads a file with peptide sequences and iRT values, and stores them into
    a dictionary irtPeptides = {"sequence1" : iRT1 , ...}
    The input file must have the following headers: sequence , iRT
    """
    irtPeptides ={}
    
    fs = csv.reader(open(filename), delimiter="\t")  
    headerset = False
    
    sequence_header = "sequence"
    irt_header      = "iRT"
    
    for cnt, srow in enumerate(fs):
        if '#' in srow[0] : continue

        #check if this row is the header & if so, set the header
        if (sequence_header in srow) and not  headerset :
            header_d = dict([ (l,i) for i,l in enumerate(srow)])
            headerset = True
            
            if irt_header not in header_d :
                print("Error: the irt peptides file does not contain a header 'iRT' !! Check the format of the irt peptides file. It should be a two columns file with the headers : sequence , iRT")
                sys.exit(2)
            continue
        
        if not headerset : continue
        
        irtPeptides[srow[header_d[sequence_header]]] = float(srow[header_d[irt_header]])
        
    return irtPeptides

def main(argv) :

    irtmodel_file   = ''
    useMinutes      = False   # [ 'seconds' , 'minutes' ]
    irtmodels       = {}
    irts_in_samples = {}
    
    useIrtPeptides  = False
    irtPeptides     = {
                       'LGGNEQVTR/2'         :   -28.3083 ,
                       'GAGSSEPVTGLDAK/2'    :   0.227424 ,
                       'VEATFGVDESNAK/2'     :   13.1078 ,
                       'YILAGVENSK/2'        :   22.3798 ,
                       'TPVISGGPYEYR/2'      :   28.9999 ,
                       'TPVITGAPYEYR/2'      :   33.6311 ,
                       'DGLDAASYYAPVR/2'     :   43.2819 ,
                       'ADVTPADFSEWSK/2'     :   54.969  ,
                       'GTFIIDPGGVIR/2'      :   71.3819 ,
                       'GTFIIDPAAVIR/2'      :   86.7152 ,
                       'LFLQFGAQGSPFLK/2'    :   98.0897
                       }
    
    #Get options
    try:
        opts, args = getopt.getopt(argv, "hi:t:ap:",["help","irtmodel","timescale","auto-align","irtpeptides"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    argsUsed = 0
    for opt,arg in opts:
        if opt in ("-h","--help") :
            usage()
            sys.exit()
        if opt in ("-a","--auto-align") :
            useIrtPeptides = True
            argsUsed += 1
        if opt in ("-i","--irt") :
            irtmodel_file = arg
            argsUsed += 2
        if opt in ("-p","--irtpeptides") :
            irtPeptides = read_irtPeptides(arg)
            argsUsed += 2
        if opt in("-t","--timescale") :
            if arg in ["minutes","seconds"] :
                if arg in ["minutes"] : useMinutes = True
            else :
                print("Choose a right time-scale. Options are: minutes, seconds")
                sys.exit(10)
            argsUsed += 2
        
    
    if len(irtmodel_file) > 0 : irtmodels = read_irtmodels(irtmodel_file, useMinutes)
    #else : 
    #    print "There is no model file!! Try again!"
        
    sptxtfiles_pat = argv[argsUsed:]
    sptxtfiles = []
    for pat in sptxtfiles_pat :
        sptxtf = glob.glob(pat)
        for file in sptxtf : sptxtfiles.append(file)


    for sptxtfile in sptxtfiles :
        print("Reading : " , sptxtfile)
        assert sptxtfile[-6:] == '.sptxt'
        if not os.path.exists(sptxtfile):
            print("The file: %s does not exist!" % sptxtfile)
            sys.exit(2)

        if useIrtPeptides :
            print("Calibration peptides :" , irtPeptides)
            irts_in_samples = search_irtPeptides(sptxtfile, irtPeptides)
            irtmodels = cal_irtmodels(irts_in_samples)

            
        #Calculate iRTs
        if len(irtmodels) > 0 : update_iRTs(sptxtfile,irtmodels)
        else :
            print("There are no models!! Review your model file.")
        
        print("done.")

if __name__ == '__main__':
    main(sys.argv[1:])
