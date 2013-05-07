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

import sys
import os
import csv
import getopt
from configobj        import ConfigObj

from msproteomicstoolslib.data_structures.modifications    import     Modifications
from msproteomicstoolslib.format.ProteinDB                import     ProteinDB  
import msproteomicstoolslib.format.speclib_db_lib        as         speclib_db_lib  


def usage() :
    print ""
    print "spectrast2tsv.py"
    print ("-" * 15)
    print "This script is used as filter from spectraST files to swath input files."
    print ""
    print "Usage: "
    print "python spectrast2tsv.py [options] spectrast_file(s)"
    print "-h                  Display this help"
    print "-d                  Remove duplicate masses from labeling"
    print "-e                  Use theoretical mass"
    print "-f    fasta_file    Fasta file to relate peptides to their proteins (this is optional)."
    print "-g    mass_modifs   List of allowed fragment mass modifications. Useful for phosphorylation and neutral losses. Example: -g -80,-98,-17,-18"
    print "-i    labeling_file File containing the amino acid isotopic labeling mass shifts. If this option is used, heavy transitions will be generated."
    print "-k    output_key    Select the output provided. Keys available: openswath, peakview. Default: peakview"
    print "-l    mass_limits   Lower and Upper mass limits. Example: -l 400,1200"
    print "-m    mods_file     File with the modifications delta mass"
    print "-n    int           Max number of reported ions per peptide/z. Default: 20"
    print "-o    int           Min number of reported ions per peptide/z. Default: 3"
    print "-p    float         Maximum error allowed at the annotation of a fragment ion. Default: 0.05"
    print "-s    ion_series    List of ion series to be used. Example: -s y,b"
    print "-t    time-scale    Options: minutes, seconds. Default: seconds."
    print "-w    swaths_file   File containing the swath ranges. This is used to remove transitions with Q3 falling in the swath mass range. (line breaks in windows/unix format)"
    print "-x    allowed_frg_z Fragment ion charge states allowed. Default: 1,2"
    print "-a    outfile       Output file name (default: appends _peakview.txt)"
    print ""

def writeStandardConfigFile(filename):

    config = ConfigObj()

    config.filename = filename

    # {'M': {'fv': ('-1.3094', '0'), 'sr': ('-13.0000', '0'), 'is': ('25.6100', '0'), 'pb': ('0.9500', '0'), 'ex': ('0.3476', '0'), 'hs': ('30.1900', '0'), 'sc': ('17.1900', '0')}}
    # {'S': {'pb': ('0.9774', '0'), 'xc': ('4.2610', '0'), 'dc': ('1.0000', '0'), 'fv': ('1.3848', '0')}}


    #Section SEQUEST
    config['SEQUEST'] = {}
    config['SEQUEST']['id'] = 'S'
    config['SEQUEST']['pb'] = 0.5
    config['SEQUEST']['xc'] = 2.5
    config['SEQUEST']['dc'] = 0
    config['SEQUEST']['fv'] = 0


    #Section MASCOT
    config['MASCOT'] = {}
    config['MASCOT']['id'] = 'M'
    config['MASCOT']['fv'] = 0
    config['MASCOT']['sr'] = 0
    config['MASCOT']['is'] = 0
    config['MASCOT']['pb'] = 0
    config['MASCOT']['ex'] = 0
    config['MASCOT']['hs'] = 0
    config['MASCOT']['sc'] = 0

    config.write()


def readLabelingFile(labeling_file) :
    #Returns a dictionary of amino-acides (including also C-Term and N-Term) with the mass shifts due to an isotope labeling experiment.
    labeling = {}

    file = open (labeling_file,"r")
    while True :
        line = file.readline()
        if len(line) == 0 : break
        if line[0] == '#' : continue

        aminoacid = ''
        mass_shift = 0.0

        sline = line.split('\t')
        if len(sline) >= 2 :
            aminoacid  = sline[0]
            mass_shift = float(sline[1])

        labeling[aminoacid] = mass_shift

    file.close()

    print "Labeling file :" , labeling
    return labeling

def read_swathsfile(swathsfile) :

    swaths = []

    fsw = open( swathsfile , 'r')

    counterline = 0

    while True:
        row = fsw.readline()

        counterline += 1
        if len(row) == 0 :
            print "Swaths file contains %s swaths" % counterline
            break

        if row[0] == '#' : continue

        srow = row.split("\t")

        if len(srow) != 2 :
            print "Error when reading swaths file. Are there more than two values in the same row?"
            sys.exit(2)

        for value in srow :
            if not is_number(value) :
                print "Error when reading swaths file. Some value(s) are not numbers!"
                sys.exit(2)

        swaths.append( ( float(srow[0]) , float(srow[1]) ) )

    return swaths

def is_Q3_in_swath_range(q1 , q3 , swaths) :
    #determine which is the swath
    swath = []
    for sw in swaths :
        if q1 >= sw[0] and q1 <= sw[1] :
            swath = sw

    #If there is no swath --> return True, so that peptides out of Q1 experiment limits are removed
    if len(swath) == 0 : return True

    if ( q3 >= swath[0] and q3<=swath[1] ) : return True
    else : return False

def removeDuplicates(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x

    seen = {}
    result = []

    for item in seq:
        marker = idfun(item)
        if marker in seen: continue

        seen[marker] = 1
        result.append(item)

    return result


def removeSimilarDuplicates(seq, tolerance , idfun=None) :
    # order preserving
    if idfun is None:
        def idfun(x): return x

    seen = {}
    result = []

    for item in seq:
        catchup = False
        marker = idfun(item)
        if not is_number(idfun(item)) : raise "comparison values are supposed to be numbers!"

        for cp in seen.iterkeys() :
            if abs( marker - cp ) <= tolerance : catchup = True

        if catchup : continue

        seen[marker] = 1
        result.append(item)

    return result

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def filterBySearchEngineParams(searchEngineInfo, parameter_thresholds) :
    spectrumOK  = True
    #print parameter_thresholds
    #print searchEngineInfo

    if not 'id' in parameter_thresholds :
        return True

    if not parameter_thresholds['id'] in searchEngineInfo :
        return True

    for parameter,threshold in parameter_thresholds.iteritems() :
        if parameter in ('id') : continue
        if parameter in searchEngineInfo[parameter_thresholds['id']] :
            if float( searchEngineInfo[parameter_thresholds['id']][parameter][0] ) < float(threshold) :
                spectrumOK = False

    return spectrumOK


def main(argv) :

    fastafile        = ''
    swathsfile        = ''
    masslimits        = [0,30000]
    ionseries        = []
    useexactmass    = False
    modificationsfile = ''
    maxtransitions    = 20
    mintransitions    = 3
    frgchargestate    = [1,2]
    precision        = 0.05
    searchEngineconfig = []
    gain_or_loss_mz_txt = []
    gain_or_loss_mz = []
    labeling = {}
    codes = ["openswath", "peakview"]
    removeDuplicatesInHeavy = False
    useMinutes = False
    keys            = ['openswath','peakview']
    key                = 'peakview'
    outputfile = None

    csv_headers_peakview =     [    'Q1', 'Q3', 'RT_detected', 'protein_name', 'isotype',
                     'relative_intensity', 'stripped_sequence', 'modification_sequence', 'prec_z',
                     'frg_type', 'frg_z', 'frg_nr', 'iRT', 'uniprot_id', 'decoy'
                     ]
    csv_headers_openswath = ['PrecursorMz', 'ProductMz', 'Tr_recalibrated', 'transition_name', 'CE',
                            'LibraryIntensity', 'transition_group_id', 'decoy', 'PeptideSequence', 'ProteinName', 
                            'Annotation', 'FullUniModPeptideName', 
                            'PrecursorCharge', 'GroupLabel', 'UniprotID', 'FragmentType', 'FragmentCharge',
                            'FragmentSeriesNumber']
    
    csv_headers = csv_headers_peakview

    swaths =[(400,425),(424,450),(449,475),(474,500),(499,525),
             (524,550),(549,575),(574,600),(599,625),(624,650),
             (649,675),(674,700),(699,725),(724,750),(749,775),
             (774,800),(799,825),(824,850),(849,875),(874,900),
             (899,925),(924,950),(949,975),(974,1000),(999,1025),
             (1024,1050),(1049,1075),(1074,1100),(1099,1125),(1124,1150),
             (1149,1175),(1174,1200)]

    #Get options
    try:
        opts, _ = getopt.getopt(argv, "hf:l:s:en:m:o:w:c:z:g:i:dx:p:t:k:a:",["help","fasta","limits","series","exact","max","modifications","min","swaths","config","writeconfig","gain","isot-labeling","remove-duplicates","charge","precision","timescale","key","output"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    argsUsed = 0
    for opt,arg in opts:
        if opt in ("-h","--help") :
            usage()
            sys.exit()
        if opt in ("-a", "--output") :
            outputfile = arg
            argsUsed += 2
        if opt in ("-f","--fasta") :
            fastafile = arg
            argsUsed += 2
        if opt in ("-m","--modifications") :
            modificationsfile = arg
            argsUsed += 2
        if opt in ("-w","--swaths") :
            print "swathsfile : " , arg
            swathsfile = arg
            argsUsed += 2
        if opt in ("-l","--limits") :
            masslimits = []
            masslimits_txt = arg.split(',')
            try :
                for val in masslimits_txt : masslimits.append( float(val) )
            except :
                print "Mass range limits are not a number! Please, try again."
                sys.exit(2)
            argsUsed += 2
        if opt in ("-s","--series") :
            ionseries = arg
            argsUsed += 2
        if opt in ("-e","--exact") :
            useexactmass = True
            argsUsed += 1
        if opt in ("-n","--max") :
            try :
                maxtransitions = int(arg)
            except :
                print "Max number of transitions is not an integer! Please, try again."
                sys.exit(2)
            argsUsed += 2
        if opt in ("-o","--min") :
            try :
                mintransitions = int(arg)
            except :
                print "Min number of transitions is not an integer! Please, try again."
            argsUsed += 2
        if opt in ("-c","--config") :
            searchEngineconfig = ConfigObj(arg)
            argsUsed+=2
        if opt in ("-z","--writeconfig") :
            writeStandardConfigFile(arg)
            sys.exit()
        if opt in ("-g","--gain") :
            gain_or_loss_mz_txt = arg.split(',')
            for obj in gain_or_loss_mz_txt :
                gain_or_loss_mz.append(float(obj))
            argsUsed += 2
        if opt in("-i","--isot-labeling") :
            labeling = readLabelingFile(arg)
            argsUsed += 2
        if opt in("-d","--remove-duplicates") :
            removeDuplicatesInHeavy = True
            argsUsed += 1
        if opt in("-x","--charge") :
            frgchargestate = []
            frgchargestate_txt = arg.split(',')
            for obj in frgchargestate_txt :
                frgchargestate.append(int(obj))
            argsUsed += 2
        if opt in("-p","--precision") :
            precision = float(arg)
            argsUsed += 2
        if opt in("-t","--timescale") :
            argsUsed += 2
            if arg in ["minutes","seconds"] :
                if arg in ["minutes"] : useMinutes = True
            else :
                print "Choose a right time-scale. Options are: minutes, seconds"
                sys.exit(10) 
        if opt in ('-k', 'key') :
            if arg not in codes :
                print "Error: key option is not valid! key : " , arg
                print "Valid options are : " , keys
                sys.exit(2)
            key = arg
            if key == 'openswath'     : csv_headers = csv_headers_openswath
            if key == 'peakview'    : csv_headers = csv_headers_peakview
            argsUsed += 2
            

    print "Masslimits:",masslimits

    if mintransitions > maxtransitions :
        print "This might seem a bit fool, but... You can't select a minimum number of transitions higher than the maximum!! "
        print "Min : " , mintransitions , " Max :" , maxtransitions
        sys.exit(2)


    sptxtfiles = argv[argsUsed:]
    if len(sptxtfiles) == 0:
        print "No input files given"
        sys.exit(2)
    
    #If a modifications file is provided, update the Modifications
    modificationsLib = Modifications()     #None
    if len(modificationsfile) > 0 :
        modificationsLib.readModificationsFile(modificationsfile)
    print "Modifications used : ",modificationsLib.mods_TPPcode.keys()

    #If a fasta file is provided, read and store it into a dictionary
    
    proteins = None
    if len(fastafile) > 0 :
        proteins = ProteinDB()
        print "Reading fasta file :" , fastafile
        proteins.readFasta(fastafile)


    #Read swaths file (if provided)
    if swathsfile != '' :
        swaths = read_swathsfile(swathsfile)
    else:
        print "Using default swath windows",swaths

    for sptxtfile in sptxtfiles :
        print "Reading : " , sptxtfile
        if not os.path.exists(sptxtfile):
            print "The file: %s does not exist!" % sptxtfile
            sys.exit(2)


        if outputfile is None:
            peakviewfilename = sptxtfile[:-6] + "_peakview.txt"
        else:
            peakviewfilename = outputfile
        try :
            writer = csv.writer(open(peakviewfilename,'w'), dialect='excel-tab')
        except :
            print "something went wrong while trying to write the file :" , peakviewfilename
            sys.exit(1)

        #write the headers
        writer.writerow( csv_headers )

        library_key = 99
        spectrastlib = speclib_db_lib.Library(library_key)

        num_spectrum = 0
        offset = spectrastlib.get_first_offset(sptxtfile)
        last_offset = -100
        modification_code = 'TPP'
        
        transition_cnt = 0
        precursor_cnt = 0
        while ( offset - last_offset > 10) :
            last_offset = offset
            offset , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)

            #for property, value in vars(spectrum).iteritems():
            #    if property    in ['compress_spectra' ] : continue
            #    print property, ": ", value
            #sys.exit()

            sequence = spectrum.name.split('/')[0]
            z_parent = float(spectrum.name.split('/')[1])

            #print sequence, z_parent

            #Declare the peptide
            pep = modificationsLib.translateModificationsFromSequence(sequence, modification_code)
            
            irt_sequence = -100
            RT_experimental = 0.0
            if spectrum.RetTime_detected != -1 :
                RT_experimental = spectrum.RetTime_detected / 60.0   #PeakView expect minutes, and spectraST reports seconds.

            if not useMinutes : RT_experimental = RT_experimental * 60


            ###only if fasta file set
            spec_proteins = []
            if proteins : spec_proteins = proteins.get_proteins_containing_peptide(pep.sequence)

            protein_code1 = ''
            protein_desc  = ''

            for prot in spec_proteins :
                protein_code1    += prot.code1
                protein_code1    += ','
                protein_desc    += prot.description
                protein_desc    += ','

            if len(protein_code1) > 0 : protein_code1 = protein_code1[:-1]
            if len(protein_desc) > 0 :  protein_desc  = protein_desc[:-1]

            if len(protein_code1) == 0 :
                if hasattr(spectrum, 'protein_ac') : protein_code1 = spectrum.protein_ac
                else : protein_code1 = 'unknown'
            if len(protein_desc)  == 0 :
                if hasattr(spectrum, 'protein_ac') : protein_desc = spectrum.protein_ac
                else : protein_desc  = 'unknown'
            ###endfasta
            
            
            num_spectrum = num_spectrum +1
            if (num_spectrum % 1000 == 0) : print "spectra processed: %s" % num_spectrum

            precursorMZ = spectrum.precursorMZ
            if useexactmass : #calculate Q1 and Q3 mass/charge values
                precursorMZ = pep.getMZ(z_parent , label = '')

            searchenginefiltered = False
            try :
                for searchengine in searchEngineconfig :
                    if not filterBySearchEngineParams(spectrum.searchEngineInfo, searchEngineconfig[searchengine] ) :
                        searchenginefiltered = True
                        continue
            except AttributeError :
                pass

            peaks = spectrum.get_peaks()
            if searchenginefiltered : peaks = []

            filteredtransitions = []
            precursor_cnt += 1
            for peak in peaks :
                if peak.is_unknown : continue
                if peak.frg_is_isotope    : continue
                if peak.frg_z not in frgchargestate : continue
                if hasattr(peak, 'mass_error') :
                    if abs(peak.mass_error) > precision : continue
                #If exact mass selected, calculate the fragment mass, otherwise keep "experimental" mass
                fragment_mz = float(peak.peak)
                if useexactmass :
                    fragment_mz  = pep.getMZfragment(peak.frg_serie , peak.frg_nr , peak.frg_z , label = '')

                #If ion series were specified by the user, filter those not matching user preferences.
                if len(ionseries) > 0 :
                    if peak.frg_serie not in ionseries : continue

                #Filter by mass range
                if fragment_mz < masslimits[0] : continue
                if fragment_mz > masslimits[1] : continue

                #Filter fragment losses and gains, isotopes...
                if peak.is_frg_loss :
                    objfound = False
                    for mz in gain_or_loss_mz :
                        #print mz, peak.frg_loss, abs(mz + peak.frg_loss[0] )
                        if abs(mz + peak.frg_loss[0]) < 0.05 :
                            objfound = True
                            if useexactmass : fragment_mz += (mz / peak.frg_z )
                            peak.frg_serie = peak.frg_serie + str(int(round(mz)))
                    if not objfound : continue
                if peak.is_frg_gain :
                    objfound = False
                    for mz in gain_or_loss_mz :
                        if abs(mz - peak.frg_gain[0]) < 0.05 :
                            objfound = True
                            if useexactmass : fragment_mz += ( mz / peak.frg_z )
                            peak.frg_serie = peak.frg_serie + '+' + str(int(round(mz)))
                    if not objfound : continue


                #Write the data into the data matrix (transitions)

                code = 'ProteinPilot'
                if key == 'openswath'     : code = 'unimod'
                if key == 'peakview'    : code = 'ProteinPilot'

                transition = []
                transition_cnt += 1
                if key == 'peakview' :
                    transition = [ precursorMZ , fragment_mz , RT_experimental , protein_desc , 'light' ,
                                    peak.intensity , spectrum.sequence , pep.getSequenceWithMods(code) , int(z_parent) ,
                                    peak.frg_serie , peak.frg_z , peak.frg_nr , irt_sequence , protein_code1 , 'FALSE']
                if key == 'openswath' :
                    transition = [precursorMZ, fragment_mz, RT_experimental, "%s_%s_%s" % (transition_cnt, pep.getSequenceWithMods(code), int(z_parent)), '-1',
                            peak.intensity, "%s_%s_%s" % (precursor_cnt, pep.getSequenceWithMods(code), int(z_parent)), 0, spectrum.sequence, protein_desc, 
                            peak.peak_annotation, pep.getSequenceWithMods(code),
                            int(z_parent), 'light', protein_code1, peak.frg_serie, peak.frg_z,
                            peak.frg_nr ]
                
                filteredtransitions.append(transition)
            
            #Sort transitions by frg_serie, frg_nr, frg_z and MINUS intensity, then remove duplicates
            #this means of every frag_serie/no/chg only the highest(!) intensity peak is stored
            filteredtransitions = sorted(filteredtransitions, key= lambda x: (x[9], x[10], x[11], -x[5]))
            filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[9], x[10], x[11]))

            #REVERSE sort the transitions by intensity
            #this means the highest(!) intensity peaks of this transition are on top
            filteredtransitions = sorted ( filteredtransitions    , key=lambda transition : transition[5], reverse=True )

            #Remove transitions with very similar Q3 masses
            massTolerance = 0.02 #This should be configured by user --> TO-DO
            filteredtransitions = removeSimilarDuplicates (filteredtransitions , massTolerance , lambda x : x[1])

            #If a swaths file was provided, remove transitions falling into the swath
            filteredtransitions_tmp = []
            if len(swaths) > 0 :
                for tr in filteredtransitions :
                    if not is_Q3_in_swath_range(tr[0] , tr[1] , swaths) : filteredtransitions_tmp.append(tr)
                filteredtransitions = filteredtransitions_tmp

            #if less transitions than the minimum --> continue to next spectrum
            #print filteredtransitions
            if len(filteredtransitions) < mintransitions :
                filteredtransitions = [] #I don't think this is really necessary, just in case.
                continue

            #Write in the peakview input file (until the max number of transitions per peptide/z (the most intense N transitions)
            for index in range(0,min(len(filteredtransitions),maxtransitions)) :
                writer.writerow(filteredtransitions[index])

                #Isotopic labeling 
                #To-Do : All the labeling calculations must be reviewed and adapted to the Peptide class
                if len(labeling) > 0 :
                    heavy_transition       = filteredtransitions[index]
                    if key == 'peakview' :
                        heavy_transition[4] = 'heavy'
                        precursorMZ_heavy      = heavy_transition[0]
                        fragment_mz_heavy      = heavy_transition[1]
                        sequence_heavy         = heavy_transition[6]
                        z_parent               = heavy_transition[8]
                        frg_serie              = heavy_transition[9]
                        frg_z                   = heavy_transition[10]
                        frg_number               = heavy_transition[11]
                    if key == 'openswath' :
                        heavy_transition[4] = 'heavy'
                        precursorMZ_heavy      = heavy_transition[0]
                        fragment_mz_heavy      = heavy_transition[1]
                        sequence_heavy         = heavy_transition[8]
                        z_parent               = heavy_transition[15]
                        frg_serie              = heavy_transition[18]
                        frg_z                   = heavy_transition[19]
                        frg_number               = heavy_transition[20]
                                        
                    #NOTE: This only works for y- and b- ions (AND their neutral losses). Other fragment series are ignored, and no heavy transition will be generated for them.
                    if frg_serie[0] not in ['y','b'] : continue


                    for aa in sequence_heavy :
                        if aa in labeling : precursorMZ_heavy += labeling[aa] / z_parent

                    frg_seq = sequence_heavy
                    #b series
                    if frg_serie == 'b': frg_seq = sequence_heavy[:frg_number]
                    #y series
                    if frg_serie == 'y': frg_seq = sequence_heavy[-frg_number:]

                    for aa in frg_seq :
                        if aa in labeling : fragment_mz_heavy += labeling[aa] / frg_z

                    #Check for C- and N-terminal labelings
                    if 'C-term' in labeling :
                        precursorMZ_heavy += labeling['C-term'] / z_parent
                        if frg_serie == 'y' : fragment_mz_heavy += labeling['C-term'] /frg_z
                    if 'N-term' in labeling :
                        precursorMZ_heavy += labeling['N-term'] / z_parent
                        if frg_serie == 'b' : fragment_mz_heavy += labeling['N-term'] /frg_z


                    #NOTE : if for any reason the Q3 mass has not changed (the labeling is not present in this fragment ion), it is not reported.
                    #        if Q3 remains the same, then the mass shift in Q1 should be enough to fall in a different SWATH window.
                    #        (In case a swaths mass ranges file is present, otherwise, we remove the not changing Q3 masses anyway).
                    if removeDuplicatesInHeavy and fragment_mz_heavy == filteredtransitions[index][1] :
                        if len(swaths) > 0 :
                            if is_Q3_in_swath_range(precursorMZ_heavy , filteredtransitions[index][0] , swaths) : continue
                        else : continue

                    #Write the heavy transition into the file
                    heavy_transition[0] = precursorMZ_heavy
                    heavy_transition[1] = fragment_mz_heavy
                    writer.writerow(heavy_transition)

            filteredtransitions = []

        print "file written : " , peakviewfilename


    print "Done."



if __name__ == '__main__':
    main(sys.argv[1:])