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
from __future__ import division
import sys
import os
import csv
import getopt
import multiprocessing
import time
import math
import sys
import numbers
from configobj import ConfigObj

from msproteomicstoolslib.data_structures.aminoacides     import Aminoacides
from msproteomicstoolslib.data_structures.modifications    import Modifications
from msproteomicstoolslib.format.ProteinDB                import ProteinDB  
import msproteomicstoolslib.format.speclib_db_lib            as speclib_db_lib  
#import msproteomicstoolslib.utils.logs.MultiProcessingLog     as MPlog 


def usage() :
    print("")
    print("spectrast2tsv.py")
    print("-" * 15)
    print("This script is used as filter from spectraST files to swath input files.")
    print("")
    print("Usage: ")
    print("python spectrast2tsv.py [options] spectrast_file(s)")
    print("-h                  Display this help")
    print("-d                  Remove duplicate masses from labeling")
    print("-e                  Use theoretical mass")
    print("-f    fasta_file    Fasta file to relate peptides to their proteins (this is optional).")
    print("-g    mass_modifs   List of allowed fragment mass modifications. Useful for phosphorylation and neutral losses. Example: -g -79.97,-97.98,-17.03,-18.01")
    print("-i    labeling_file File containing the amino acid isotopic labeling mass shifts. If this option is used, heavy transitions will be generated.")
    print("-k    output_key    Select the output provided. Keys available: openswath, peakview. Default: peakview")
    print("-l    mass_limits   Lower and upper mass limits of fragment ions. Example: -l 400,2000")
    print("-m    mods_file     File with the modifications delta mass")
    print("-n    int           Max number of reported ions per peptide/z. Default: 20")
    print("-o    int           Min number of reported ions per peptide/z. Default: 3")
    print("-p    float         Maximum error allowed at the annotation of a fragment ion. Default: 0.05")
    print("-q    int            Number of processors to use (only for isoforms!). Default: 1")
    print("-s    ion_series    List of ion series to be used. Example: -s y,b")
    print("-t    time-scale    Options: minutes, seconds. Default: seconds.")
    print("-u     unimod-code    Use this unimod code as a switching modification. Useful for phosphorylations. Example: -u 21")
    print("-v                  Verbose mode.")
    print("-w    swaths_file   File containing the swath ranges. This is used to remove transitions with Q3 falling in the swath mass range. (line breaks in windows/unix format)")
    print("-x    allowed_frg_z Fragment ion charge states allowed. Default: 1,2")
    print("-y    UIS-order     When using a switching modification, this determines the UIS order to be calculated. If -1 is set, all transitions for each isoform will be reported. Default : 2")
    print("-a    outfile       Output file name (default: appends _peakview.txt)")
    print("")

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

class Label(object):

    def __init__(self, AA, deltamass, name):
        self.AA = AA
        self.deltamass = deltamass
        self.name = name

def readLabelingFile(labeling_file) :
    #Returns a dictionary of amino-acides (including also C-Term and N-Term) with the mass shifts due to an isotope labeling experiment.
    labeling = {}

    cur_file = open (labeling_file,"r")
    for line in cur_file:
        # Skip empty lines or comments
        if len(line.strip()) == 0 : continue
        if line[0] == '#' : continue

        labelname = ''
        aminoacid = ''
        mass_shift = 0.0

        sline = line.strip().split('\t')
        if len(sline) >= 2 :
            aminoacid  = sline[0]
            mass_shift = float(sline[1])
            if len(sline) >= 3:
                labelname = sline[2]

        labeling[aminoacid] = Label(aminoacid, mass_shift, labelname)

    cur_file.close()

    print("Labeling file :" , labeling)
    return labeling

def read_swathsfile(swathsfile) :
    swaths = []
    with open(swathsfile) as infile:
        for count, row in enumerate(infile):
            if not row:
                print("Swaths file contains %s swaths" % count)
                break

            if row[0] == '#':
                continue

            srow = row.split("\t")
            if len(srow) != 2:
                print("Error when reading swaths file. Are there more than two values in the same row?")
                sys.exit(2)

            try:
                swaths.append((float(srow[0]), float(srow[1])))
            except ValueError as e:
                print("Error while reading swaths file: Some value(s) "
                      "are not numbers: " + str(e))
                sys.exit(2)

    return swaths
    

def is_Q3_in_swath_range(q1 , q3 , swaths):
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


def removeSimilarDuplicates(seq, tolerance, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x

    seen = set()
    result = []

    for item in seq:
        marker = idfun(item)
        if not isinstance(marker, numbers.Real):
            raise ValueError("comparison values are supposed to be numbers!")

        if any(abs(marker - cp) <= tolerance for cp in seen):
            continue

        seen.add(marker)
        result.append(item)

    return result


def filterBySearchEngineParams(searchEngineInfo, parameter_thresholds) :
    spectrumOK  = True
    #print parameter_thresholds
    #print searchEngineInfo

    if not 'id' in parameter_thresholds :
        return True

    if not parameter_thresholds['id'] in searchEngineInfo :
        return True

    for parameter,threshold in parameter_thresholds.items() :
        if parameter in ('id') : continue
        if parameter in searchEngineInfo[parameter_thresholds['id']] :
            if float( searchEngineInfo[parameter_thresholds['id']][parameter][0] ) < float(threshold) :
                spectrumOK = False

    return spectrumOK


def get_iso_species(sptxtfile, switchingModification, modificationsLib, aaLib = None) :
    '''It reads a spectraST file, and gets a dictionary of peptides containing a dictionary of isoforms of the peptide. 
    The values of the latter dictionary are the file offset in which the related species are located at the spectraST file. If
    the specie is not located (it has only been theoretically defined), the offset value is -1. 
    Example : iso_species = { 'PEPTSISE_(UniMod:21)' : {'PEPTS(UniMod:21)ISE' : 335542366 , 'PEPTSIS(UniMod:21)E' : -1 } },
    where PEPTSISE_(UniMod:21) means that the peptide contains only one phosphorilation. Two phosphorilations are reported as 
    PEPTSISE_(UniMod:21)(UniMod:21)
    '''
    iso_species = {}
    library_key = 99
    spectrastlib = speclib_db_lib.Library(library_key)
    num_spectrum = 0
    offset = spectrastlib.get_first_offset(sptxtfile)
    last_offset = -100
    unimodcode = ''
    switchingMod = modificationsLib.mods_unimods[switchingModification]
    unimodcode = switchingMod.getcode('unimod')[1:]

    while ( offset - last_offset > 10) :
        last_offset = offset
        offset , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)
        sequence = spectrum.name.split('/')[0]
        z_parent = float(spectrum.name.split('/')[1])
        #Get the possible modification switching sites. If there is only one possibility, discard the peptide for the dictionary
        pep = modificationsLib.translateModificationsFromSequence(sequence, 'TPP', aaLib = aaLib)
        isobaric_peptides = []
        if len(pep.modifications) == 0 : continue
        #Count the number of switching mods
        num_of_swmods = 0
        for pos,mod in pep.modifications.items() :
            if mod.unimodAccession == switchingModification : num_of_swmods += 1
        #Count the number of modification sites in the peptide
        aaTargets = [ m.aminoacid for m in modificationsLib.list if m.unimodAccession == switchingModification ]
        aaList = pep._getAminoacidList()
        num_of_sites = 0
        for aa in aaTargets :
            if aa in aaList : num_of_sites += aaList[aa]
        if num_of_sites < 2 : continue 
        #If the peptide+#switchingMods is not present in the dictionary, initialize its key with all the possible isobaric species.
        #Otherwise, find the corresponding specie in the dictionary, and change its offset value.
        specie_family = pep.sequence + "_%s" % (num_of_swmods * unimodcode)
        specie = pep.getSequenceWithMods('unimod')
        if specie_family in iso_species : 
            if specie in iso_species[specie_family] : iso_species[specie_family][specie] = last_offset
        else :
            iso_species[specie_family] = {}
            #Calculate here all the possible isobaric peptides, and initialize the value in the dictionary
            isobaric_peptides = pep.calIsoforms(switchingMod,modificationsLib)    
            for isopep in isobaric_peptides : 
                iso_species[specie_family][isopep.getSequenceWithMods('unimod')] = -1
                if isopep.getSequenceWithMods('unimod') == specie : iso_species[specie_family][specie] = last_offset 
    return iso_species

def isoform_writer(isobaric_species, lock, sptxtfile, modLibrary, aaLib, searchEngineconfig, masslimits, switchingModification_, 
                writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key, useMinutes, precision):
    filteredtransitions = []
    library_key = 99
    spectrastlib = speclib_db_lib.Library(library_key)
    switchingModification = modLibrary.mods_unimods[switchingModification_]
    
    pepfamily_cnt = 0
    precursor_cnt = 0   # To-Do : This maybe I should pass to the worker the indexes, so that we don't repeat indexes. 
    transition_cnt = 0  # To-Do : This maybe I should pass to the worker the indexes, so that we don't repeat indexes.
    fut_progress = 0.05
    for pepfamily , isoforms in isobaric_species.items() :
        pepfamily_cnt += 1
        progress = pepfamily_cnt / len(isobaric_species)
        if progress >= fut_progress and progress > 0.005 : 
            print("process id:", os.getpid() , ", ", fut_progress*100 ,  "% of peptide families written.")
            fut_progress += 0.05 
        
        #print pepfamily, isoforms
        precursor_cnt  += 1
        protein_code1 = pepfamily
        protein_desc  = pepfamily
        #for isoform, offset in isoforms.iteritems() :
        #    print "%s : %s"  % (isoform , offset)
        for isoform_,offset in isoforms.items() :
            spectrum = None
            isoform  = modLibrary.translateModificationsFromSequence(isoform_, 'unimod', aaLib = aaLib)
            
            if offset > 0 :  
                if lock : lock.acquire()
                _ , spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)
                time.sleep(0.001)
                if lock : lock.release()
            #Get the shared and unshared ions of this isoform to all its family
            otherIsoforms = []
            for isof in isoforms:
                if isof != isoform_:
                    otherIsoforms.append(modLibrary.translateModificationsFromSequence(isof, 'unimod', aaLib=aaLib))
            
            shared, unshared = isoform.comparePeptideFragments(otherIsoforms, ['y','b'], precision = 1e-5)
        
            #if there is a spectrum for the isoform, use it. Otherwise, use dummy data
            z_parent = 2  #To-Do: Better if we'd take the most common parental charge among the family
            sequence = isoform.getSequenceWithMods('TPP') # spectrum.name.split('/')[0]
            iRT_experimental = -100000.0
            RT_experimental = -100000.0
            searchenginefiltered = False
            peaks = []
            if spectrum : 
                z_parent = float(spectrum.name.split('/')[1])
                if spectrum.RetTime_detected:
                    RT_experimental = spectrum.RetTime / 60.0   #PeakView expect minutes, and spectraST reports seconds.
                if spectrum.iRT_detected:
                    iRT_experimental = spectrum.iRT / 60.0   #PeakView expect minutes, and spectraST reports seconds.
                else:
                    iRT_experimental = RT_experimental
                if not useMinutes : iRT_experimental = iRT_experimental * 60
                if not useMinutes : RT_experimental = RT_experimental * 60
                try :
                    for searchengine in searchEngineconfig :
                        if not filterBySearchEngineParams(spectrum.searchEngineInfo, searchEngineconfig[searchengine] ) :
                            searchenginefiltered = True
                            continue
                except AttributeError :
                    pass
                peaks = spectrum.get_peaks()
            precursorMZ = isoform.getMZ(z_parent, label ='')
            

            if searchenginefiltered : peaks = []
            
            
            for (frg_serie,frg_nr,frg_z,gainloss,fragment_mz) in unshared :
                #Check whether we have the fragment in the spectrum    
                rel_intensity = 1 #Dummy value in case there's no spectrum
                is_in_spectrum = False
                for peak in peaks :
                    ''' These filters are not useful in this use case, since they will be filtered in a previous step
                    if peak.is_unknown : continue
                    if peak.frg_is_isotope    : continue
                    if peak.frg_z not in frgchargestate : continue
                    if hasattr(peak, 'mass_error') :
                        if abs(peak.mass_error) > precision : continue
                    '''
                    if abs(float(peak.peak) - fragment_mz) < precision : 
                        is_in_spectrum = True 
                        rel_intensity = float(peak.intensity)
                        break
                #print isoform.getSequenceWithMods('unimod') , frg_serie, frg_nr, frg_z, fragment_mz, is_in_spectrum, rel_intensity
                #Filter by mass range
                if fragment_mz < masslimits[0] : continue
                if fragment_mz > masslimits[1] : continue
                    
                code = 'ProteinPilot'
                if key == 'openswath'     : code = 'unimod'
                if key == 'peakview'     : code = 'ProteinPilot'

                transition = []
                transition_cnt += 1
                if key == 'peakview' :
                    transition = [ precursorMZ , fragment_mz , iRT_experimental , protein_desc , 'light' ,
                                    rel_intensity , isoform.sequence , isoform.getSequenceWithMods(code) , int(z_parent) ,
                                    frg_serie , frg_z , frg_nr , iRT_experimental , protein_code1 , 'FALSE']
                if key == 'openswath' :
                    transition = [precursorMZ, fragment_mz, iRT_experimental, "%s_%s_%s" % (transition_cnt, isoform.getSequenceWithMods(code), int(z_parent)), '-1',
                            rel_intensity, "%s_%s_%s" % (precursor_cnt, isoform.getSequenceWithMods(code), int(z_parent)), 0, isoform.sequence, protein_desc, 
                            "%s_%s_%s" %(frg_serie, frg_z, frg_nr), isoform.getSequenceWithMods(code), int(z_parent), 'light', protein_code1, frg_serie, frg_z, frg_nr ]
                filteredtransitions.append(transition)
            
            #For each isoform, filtering and write
            do_filtering_and_write(filteredtransitions, writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, 0.02, lock = lock)
            filteredtransitions = []
            


def mp_isoform_writer(isobaric_species, sptxtfile, modLibrary, aaLib, searchEngineconfig, masslimits, switchingModification_, 
                writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key,useMinutes, precision, nprocs):
    def worker(pepfamilies, lock, sptxtfile, modLibrary, aaLib, searchEngineconfig, masslimits, switchingModification_, 
                writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key, useMinutes, precision):
        """ The worker function, invoked in a process. 'nums' is a
            list of pep families to process.
        """
        #outdict = {}
        #To-Do: This worker is absolutely absurd --> REDEFINE!!!!
        isoform_writer(pepfamilies, lock, sptxtfile, modLibrary, aaLib, searchEngineconfig, masslimits, switchingModification_, 
                writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key, useMinutes, precision)
        #out_q.put(outdict)

    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    #out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(len(isobaric_species) / float(nprocs)))
    print("%s peptide families divided into %s chunks of %s" % (len(isobaric_species) , nprocs, chunksize) )
    procs = []
    lock = multiprocessing.Lock()
    #lock = None
    
    dict_chunks = {}
    last_chunk = 0
    dict_chunks[last_chunk] = {}
    for idx, (pepf, value) in enumerate(isobaric_species.items()) : #divide the dictionary into nprocs different dictionaries
        if idx // chunksize > last_chunk:
            last_chunk += 1
            dict_chunks[last_chunk] = {}
        dict_chunks[last_chunk][pepf] = value

    for i in dict_chunks.values():
        p = multiprocessing.Process(
                target=worker,
                args=(i, lock, sptxtfile, modLibrary, aaLib, 
                    searchEngineconfig, masslimits, switchingModification_, 
                writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key, useMinutes, precision))
        procs.append(p)
        p.start()
        #out_q.put(p)

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    #resultdict = {}
    #for i in range(nprocs):
    #    resultdict.update(out_q.get())

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return 0


def transitions_isobaric_peptides(isobaric_species , sptxtfile, switchingModification_, modLibrary,UISorder, ionseries, fragmentlossgains, frg_z_list,  
                        useMinutes, searchEngineconfig, masslimits, key, precision,  
                        writer,  labeling, removeDuplicatesInHeavy, swaths, mintransitions, maxtransitions, massTolerance, aaLib = None, nprocs = 0, verbose = False) :
    '''
    This returns the (specific!) transitions of all the isoforms of the peptides present in the spectraST library. If there is no
    spectraST reference for an isoform, transitions are chosen randomly, respecting the established filter rules, and a dummy 
    relative intensity value is given.
    '''
    library_key = 99
    spectrastlib = speclib_db_lib.Library(library_key)
    switchingModification = modLibrary.mods_unimods[switchingModification_]
    
    filteredtransitions = []
    transition_cnt = 0
    precursor_cnt  = 0
    
    print("A total of %s peptide families have been found."  % len(isobaric_species) )
    
    pepfamily_cnt = 0
    fut_progress = 0.01

    #Multiprocess this part of the code (only if multiprocess option has been specified)
    if nprocs > 0 :
        mp_isoform_writer(isobaric_species, sptxtfile, modLibrary, aaLib, searchEngineconfig, masslimits, switchingModification_, 
                    writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, key, useMinutes, precision, nprocs)
        print("done!")
        sys.exit()
    
    for pepfamily , isoforms in isobaric_species.items() :
        #print pepfamily, [ isof for isof in isoforms]
        pepfamily_cnt += 1
        progress = pepfamily_cnt / len(isobaric_species)
        if progress >= fut_progress and progress > 0.005 : 
            print(fut_progress*100 ,  "% of peptide families written.")
            fut_progress += 0.01 
        
        #print pepfamily, isoforms
        precursor_cnt  += 1
        protein_code1 = pepfamily
        protein_desc  = pepfamily
        #for isoform, offset in isoforms.iteritems() :
        #    print "%s : %s"  % (isoform , offset)
        for isoform_,offset in isoforms.items() :
            spectrum = None
            isoform  = modLibrary.translateModificationsFromSequence(isoform_, 'unimod', aaLib = aaLib)
            
            if offset > 0:
                _, spectrum = spectrastlib.read_sptxt_with_offset(sptxtfile,offset)
            #Get the shared and unshared ions of this isoform to all its family
            otherIsoforms = []
            for isof in isoforms:
                if isof != isoform_:
                    otherIsoforms.append(modLibrary.translateModificationsFromSequence(isof, 'unimod', aaLib = aaLib))
            
            uis_list , uis_annotated_list = isoform.cal_UIS(otherIsoforms, UISorder = UISorder,  ionseries = ionseries, 
                                            fragmentlossgains = fragmentlossgains, precision = massTolerance, frg_z_list = frg_z_list, mass_limits = masslimits)
            
            #if there is a spectrum for the isoform, use it. Otherwise, use dummy data
            z_parent = 2  #To-Do: Better if we'd take the most common parental charge among the family
            sequence = isoform.getSequenceWithMods('TPP') # spectrum.name.split('/')[0]
            iRT_experimental = -100000.0
            RT_experimental = -100000.0
            searchenginefiltered = False
            peaks = []
            if spectrum : 
                z_parent = float(spectrum.name.split('/')[1])
                if spectrum.RetTime_detected:
                    RT_experimental = spectrum.RetTime / 60.0   #PeakView expect minutes, and spectraST reports seconds.
                if spectrum.iRT_detected:
                    iRT_experimental = spectrum.iRT / 60.0   #PeakView expect minutes, and spectraST reports seconds.
                else:
                    iRT_experimental = RT_experimental
                if not useMinutes : iRT_experimental = iRT_experimental * 60
                if not useMinutes : RT_experimental = RT_experimental * 60
                try :
                    for searchengine in searchEngineconfig :
                        if not filterBySearchEngineParams(spectrum.searchEngineInfo, searchEngineconfig[searchengine] ) :
                            searchenginefiltered = True
                            continue
                except AttributeError :
                    pass
                peaks = spectrum.get_peaks()
            precursorMZ = isoform.getMZ(z_parent, label ='')
            

            if searchenginefiltered : peaks = []
            #print uis_annotated_list
            for uis_masses, uis in uis_annotated_list.items() :
                for uis_unit in uis :
                    # uis_annotated_list = [('y', 6, 1, 0, 646.3406318939999), ('y', 8, 1, 0, 928.365933736)]
                    #print uis_unit
                    uis_order   = len(uis)
                    uis_serie   = uis_unit[0]
                    uis_nr      = uis_unit[1]
                    uis_frg_z   = uis_unit[2]
                    uis_lossgain= uis_unit[3]
                    uis_mass    = uis_unit[4]
                    rel_intensity = 1 #Dummy value in case there's no spectrum
                    is_in_spectrum = False
                    for peak in peaks : 
                        if abs(float(peak.peak) - uis_mass) < precision:
                            is_in_spectrum = True 
                            rel_intensity = float(peak.intensity)
                            break
                        
                    
                    if uis_mass < masslimits[0]:
                        continue   # I don't think this is necessary,
                    if uis_mass > masslimits[1]:
                        continue   # since UIS are already mass limited, but...
               
                    code = 'ProteinPilot'
                    if key == 'openswath':
                        code = 'unimod'
                    if key == 'peakview':
                        code = 'ProteinPilot'
    
                    transition = []
                    transition_cnt += 1
                    if key == 'peakview':
                        if abs(uis_lossgain) > 0.05:
                            uis_serie = uis_serie + str(int(round(uis_lossgain)))
                        transition = [
                            precursorMZ, uis_mass, iRT_experimental,
                            protein_desc, 'light', rel_intensity,
                            isoform.sequence, isoform.getSequenceWithMods(code),
                            int(z_parent), uis_serie, uis_frg_z, uis_nr,
                            iRT_experimental, protein_code1, 'FALSE',
                            uis_order, uis_masses
                        ]
                    if key == 'openswath' :
                        transition = [precursorMZ, uis_mass, iRT_experimental, "%s_%s_%s" % (transition_cnt, isoform.getSequenceWithMods(code), int(z_parent)), '-1',
                                rel_intensity, "%s_%s_%s" % (precursor_cnt, isoform.getSequenceWithMods(code), int(z_parent)), 0, isoform.sequence, protein_desc, 
                                "%s_%s_%s" %(uis_serie, uis_frg_z, uis_nr), isoform.getSequenceWithMods(code), int(z_parent), 'light', protein_code1, uis_serie, uis_frg_z, uis_nr, uis_order, uis_masses ]
                    filteredtransitions.append(transition)
            
            #For each isoform, filtering and write
            if verbose and ( len(filteredtransitions) > maxtransitions or len(filteredtransitions) ) < mintransitions : print(mintransitions, maxtransitions, len(filteredtransitions))
            do_filtering_and_write(filteredtransitions, writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, 0.02, verbose = verbose)
            filteredtransitions = []

    #print  filteredtransitions
    return 0

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
    switchingModification = None
    aaLib = Aminoacides()
    nprocs = 0
    verbose = False
    UISorder = 2
    
    csv_headers_peakview =     [    'Q1', 'Q3', 'RT_detected', 'protein_name', 'isotype',
                     'relative_intensity', 'stripped_sequence', 'modification_sequence', 'prec_z',
                     'frg_type', 'frg_z', 'frg_nr', 'iRT', 'uniprot_id', 'decoy' , 'N', 'confidence', 'shared'
                     ]
    csv_headers_openswath = ['PrecursorMz', 'ProductMz', 'Tr_recalibrated', 'transition_name', 'CE',
                            'LibraryIntensity', 'transition_group_id', 'decoy', 'PeptideSequence', 'ProteinName', 
                            'Annotation', 'FullUniModPeptideName', 
                            'PrecursorCharge', 'PeptideGroupLabel',   'UniprotID', 'FragmentType', 'FragmentCharge',
                            'FragmentSeriesNumber', 'LabelType']
    
    csv_headers = csv_headers_peakview

    swaths = []
    #swaths =[(400,425),(424,450),(449,475),(474,500),(499,525),
    #         (524,550),(549,575),(574,600),(599,625),(624,650),
    #         (649,675),(674,700),(699,725),(724,750),(749,775),
    #         (774,800),(799,825),(824,850),(849,875),(874,900),
    #         (899,925),(924,950),(949,975),(974,1000),(999,1025),
    #         (1024,1050),(1049,1075),(1074,1100),(1099,1125),(1124,1150),
    #         (1149,1175),(1174,1200)]

    #Get options
    try:
        opts, _ = getopt.getopt(argv, "hf:l:s:en:m:o:w:c:z:g:i:dx:p:t:k:a:u:q:vy:",["help","fasta","limits","series","exact","max","modifications","min","swaths","config","writeconfig","gain","isot-labeling","remove-duplicates","charge","precision","timescale","key","output","switchingmod","nprocs","verbose","UISorder"])

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
            print("swathsfile : " , arg)
            swathsfile = arg
            argsUsed += 2
        if opt in ("-l","--limits") :
            masslimits = []
            masslimits_txt = arg.split(',')
            try :
                for val in masslimits_txt : masslimits.append( float(val) )
            except :
                print("Mass range limits are not a number! Please, try again.")
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
                print("Max number of transitions is not an integer! Please, try again.")
                sys.exit(2)
            argsUsed += 2
        if opt in ("-o","--min") :
            try :
                mintransitions = int(arg)
            except :
                print("Min number of transitions is not an integer! Please, try again.")
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
                print("Choose a right time-scale. Options are: minutes, seconds")
                sys.exit(10) 
        if opt in ('-k', 'key') :
            if arg not in codes :
                print("Error: key option is not valid! key : " , arg)
                print("Valid options are : " , keys)
                sys.exit(2)
            key = arg
            if key == 'openswath'     : csv_headers = csv_headers_openswath
            if key == 'peakview'    : csv_headers = csv_headers_peakview
            argsUsed += 2
        if opt in ('-u', 'switchingmod') :
            argsUsed += 2
            switchingModification = int(arg)
        if opt in ('-q', '--nprocs') :
            argsUsed += 2
            nprocs = int(arg)
        if opt in ('-v','--verbose') :
            argsUsed += 1
            verbose = True
        if opt in ('-y','--UISorder') :
            argsUsed += 2
            UISorder = int(arg)
            

    print("Masslimits:",masslimits)

    if mintransitions > maxtransitions :
        print("This might seem a bit fool, but... You can't select a minimum number of transitions higher than the maximum!! ")
        print("Min : " , mintransitions , " Max :" , maxtransitions)
        sys.exit(2)


    sptxtfiles = argv[argsUsed:]
    if len(sptxtfiles) == 0:
        print("No input files given")
        sys.exit(2)
    
    #If a modifications file is provided, update the Modifications
    modificationsLib = Modifications()     #None
    if len(modificationsfile) > 0 :
        modificationsLib.readModificationsFile(modificationsfile)
    print("Modifications used : ", list(modificationsLib.mods_TPPcode.keys()))
    
    
    #If a fasta file is provided, read and store it into a dictionary
    
    proteins = None
    if len(fastafile) > 0 :
        proteins = ProteinDB()
        print("Reading fasta file :" , fastafile)
        proteins.readFasta(fastafile)

    protein_cnt = 1
    protein_index = {}

    #Read swaths file (if provided)
    if swathsfile != '' :
        swaths = read_swathsfile(swathsfile)
    else:
        #print "Using default swath windows",swaths
        print("No swath windows set.")

    for sptxtfile in sptxtfiles :
        print("Reading : " , sptxtfile)
        if not os.path.exists(sptxtfile):
            print("The file: %s does not exist!" % sptxtfile)
            sys.exit(2)

        if outputfile is None:
            peakviewfilename = sptxtfile[:-6] + "_peakview.txt"
        else:
            peakviewfilename = outputfile
        try :
            writer = csv.writer(open(peakviewfilename,'w'), dialect='excel-tab')
        except :
            print("something went wrong while trying to write the file :" , peakviewfilename)
            sys.exit(1)

        #write the headers
        if switchingModification : 
            csv_headers.extend(["UIS_order", "UIS_mass_list"])
        writer.writerow( csv_headers )

        #If a switching modification has been defined, a dictionary containing isobaric species must be created
        if switchingModification :
            #Check whether the switching modification is actually present in the modifications lib.
            if switchingModification not in modificationsLib.mods_unimods :
                raise Exception("Error: the switching modification given by the user is not in the modifications library!")
            switchingMod = modificationsLib.mods_unimods[switchingModification]
            isobaric_species = get_iso_species(sptxtfile, switchingModification, modificationsLib, aaLib = aaLib)
            if 0 not in gain_or_loss_mz : gain_or_loss_mz.append(0.0)
            transitions_isobaric_peptides(isobaric_species , sptxtfile, switchingModification, modificationsLib, UISorder, ionseries, gain_or_loss_mz, frgchargestate,
                        useMinutes, searchEngineconfig, masslimits, key, precision, writer,  labeling, removeDuplicatesInHeavy, swaths,mintransitions, maxtransitions, 0.02, aaLib = aaLib, nprocs = nprocs, verbose = verbose)
            
            print("done!")
            sys.exit()

        filteredtransitions = []

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
            pep = modificationsLib.translateModificationsFromSequence(sequence, modification_code, aaLib = aaLib)
            
            iRT_experimental = 0.0
            RT_experimental = 0.0
            if spectrum.RetTime_detected:
                RT_experimental = spectrum.RetTime / 60.0   #PeakView expect minutes, and spectraST reports seconds.
            if spectrum.iRT_detected:
                iRT_experimental = spectrum.iRT / 60.0   #PeakView expect minutes, and spectraST reports seconds.
            else:
                iRT_experimental = RT_experimental

            if not useMinutes : iRT_experimental = iRT_experimental * 60
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
            if (num_spectrum % 1000 == 0) : print("spectra processed: %s" % num_spectrum)

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
                            if useexactmass : fragment_mz += mz / peak.frg_z
                            peak.frg_serie = peak.frg_serie + str(int(round(mz)))
                    if not objfound : continue
                if peak.is_frg_gain :
                    objfound = False
                    for mz in gain_or_loss_mz :
                        if abs(mz - peak.frg_gain[0]) < 0.05 :
                            objfound = True
                            if useexactmass : fragment_mz += mz / peak.frg_z
                            peak.frg_serie = peak.frg_serie + '+' + str(int(round(mz)))
                    if not objfound : continue


                #Write the data into the data matrix (transitions)

                code = 'ProteinPilot'
                if key == 'openswath'     : code = 'unimod'
                if key == 'peakview'    : code = 'ProteinPilot'

                if protein_desc not in protein_index:
                	protein_index[protein_desc] = protein_cnt
                	protein_cnt+=1

                if protein_desc[:2] == "1/":
                	protein_shared = "FALSE"
                else:
                	protein_shared = "TRUE"

                transition = []
                transition_cnt += 1
                if key == 'peakview' :
                    transition = [ precursorMZ , fragment_mz , iRT_experimental , protein_desc , 'light' ,
                                    peak.intensity , spectrum.sequence , pep.getSequenceWithMods(code) , int(z_parent) ,
                                    peak.frg_serie , peak.frg_z , peak.frg_nr , iRT_experimental , protein_code1 , 'FALSE', protein_index[protein_desc] , 1 , protein_shared ]
                if key == 'openswath' :
                    transition_group_id = "%s_%s_%s" % (precursor_cnt, pep.getSequenceWithMods(code), int(z_parent))
                    transition = [precursorMZ, fragment_mz, iRT_experimental, "%s_%s%s_%s_%s_%s" % (transition_cnt, peak.frg_serie, peak.frg_nr, peak.frg_z, pep.getSequenceWithMods(code), int(z_parent)), '-1',
                            peak.intensity, transition_group_id, 0, spectrum.sequence, protein_desc, 
                            peak.peak_annotation, pep.getSequenceWithMods(code),
                            int(z_parent), transition_group_id, protein_code1, peak.frg_serie, peak.frg_z,
                            peak.frg_nr , 'light']
                
                filteredtransitions.append(transition)
            

            #Sort transitions by frg_serie, frg_nr, frg_z and MINUS intensity, then remove duplicates
            #this means of every frag_serie/no/chg only the highest(!) intensity peak is stored
            if key == 'peakview':
                filteredtransitions = sorted(filteredtransitions, key=lambda x: (x[9], x[10], x[11], -x[5]))
                filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[9], x[10], x[11]))
            elif key == 'openswath':
                filteredtransitions = sorted(filteredtransitions, key=lambda x: (x[15], x[16], x[17], -x[5]))
                filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[15], x[16], x[17]))

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
                        # heavy_transition[13] = 'heavy'
                        heavy_transition[18] = 'heavy'
                        heavy_transition[13] = heavy_transition[6]
                        precursorMZ_heavy      = heavy_transition[0]
                        fragment_mz_heavy      = heavy_transition[1]
                        sequence_heavy         = heavy_transition[8]
                        z_parent               = heavy_transition[12]
                        frg_serie              = heavy_transition[15]
                        frg_z                   = heavy_transition[16]
                        frg_number               = heavy_transition[17]

                        # Modify full sequence with name of label
                        fullseq = heavy_transition[11]
                        fullseq_n = ""
                        for aa in fullseq:
                            fullseq_n += aa
                            if aa in labeling and labeling[aa].name:
                                fullseq_n += "(%s)" % labeling[aa].name

                        # Store modified sequence and make transition id unique again
                        heavy_transition[11] = fullseq_n
                        heavy_transition[3] += 'heavy'
                        heavy_transition[6] += 'heavy'




                                        
                    # TODO also change full unimod peptide name     


                    #NOTE: This only works for y- and b- ions (AND their neutral losses). Other fragment series are ignored, and no heavy transition will be generated for them.
                    if frg_serie[0] not in ['y','b'] : continue


                    for aa in sequence_heavy :
                        if aa in labeling : precursorMZ_heavy += labeling[aa].deltamass / z_parent

                    frg_seq = sequence_heavy
                    #b series
                    if frg_serie == 'b': frg_seq = sequence_heavy[:frg_number]
                    #y series
                    if frg_serie == 'y': frg_seq = sequence_heavy[-frg_number:]

                    for aa in frg_seq :
                        if aa in labeling : fragment_mz_heavy += labeling[aa].deltamass / frg_z

                    #Check for C- and N-terminal labelings
                    if 'C-term' in labeling:
                        precursorMZ_heavy += labeling['C-term'].deltamass / z_parent
                        if frg_serie == 'y':
                            fragment_mz_heavy += labeling['C-term'].deltamass / frg_z
                    if 'N-term' in labeling:
                        precursorMZ_heavy += labeling['N-term'].deltamass / z_parent
                        if frg_serie == 'b':
                            fragment_mz_heavy += labeling['N-term'].deltamass / frg_z


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

        print("file written : " , peakviewfilename)


    print("Done.")






def do_filtering_and_write(filteredtransitions, writer, labeling, removeDuplicatesInHeavy, swaths, mintransitions, maxtransitions, massTolerance = 0.02, verbose = False , lock = None) :
    #Sort transitions by frg_serie, frg_nr, frg_z and MINUS intensity, then remove duplicates
    #this means of every frag_serie/no/chg only the highest(!) intensity peak is stored
    if verbose : print("initial number transitions " , len(filteredtransitions))
    if key == 'peakview':
        filteredtransitions = sorted(filteredtransitions, key=lambda x: (x[9], x[10], x[11], -x[5]))
        filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[9], x[10], x[11]))
    elif key == 'openswath':
        filteredtransitions = sorted(filteredtransitions, key=lambda x: (x[15], x[16], x[17], -x[5]))
        filteredtransitions = removeDuplicates(filteredtransitions, lambda x: (x[15], x[16], x[17]))
    if verbose : print("after removing duplicates " , len(filteredtransitions))

    #REVERSE sort the transitions by intensity
    #this means the highest(!) intensity peaks of this transition are on top
    filteredtransitions = sorted ( filteredtransitions    , key=lambda transition : transition[5], reverse=True )

    #Remove transitions with very similar Q3 masses
    massTolerance = 0.02 #This should be configured by user --> TO-DO
    filteredtransitions = removeSimilarDuplicates (filteredtransitions , massTolerance , lambda x : x[1])
    if verbose : print("after removing similar duplicates " , len(filteredtransitions))

    #If a swaths file was provided, remove transitions falling into the swath
    filteredtransitions_tmp = []
    if len(swaths) > 0 :
        for tr in filteredtransitions :
            if not is_Q3_in_swath_range(tr[0] , tr[1] , swaths) : filteredtransitions_tmp.append(tr)
        filteredtransitions = filteredtransitions_tmp
    if verbose : print("after removing transitions within the swath isolation window " , len(filteredtransitions))

    #if less transitions than the minimum --> continue to next spectrum
    #print filteredtransitions
    if len(filteredtransitions) < mintransitions :
        filteredtransitions = [] #I don't think this is really necessary, just in case.
    if verbose : print("after removing num of transitions below the required minimum " , len(filteredtransitions))

    adssfsfadsf
    
    #Write in the peakview input file (until the max number of transitions per peptide/z (the most intense N transitions)
    for index in range(0,min(len(filteredtransitions),maxtransitions)) :
        if lock : lock.acquire()
        writer.writerow(filteredtransitions[index])
        time.sleep(0.001)
        if lock : lock.release()
        
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
                print(heavy_transition)
                heavy_transition[4] = 'heavy'
                print(heavy_transition)
                precursorMZ_heavy      = heavy_transition[0]
                fragment_mz_heavy      = heavy_transition[1]
                sequence_heavy         = heavy_transition[8]
                z_parent               = heavy_transition[12]
                frg_serie              = heavy_transition[15]
                frg_z                   = heavy_transition[16]
                frg_number               = heavy_transition[17]
                                
            #NOTE: This only works for y- and b- ions (AND their neutral losses). Other fragment series are ignored, and no heavy transition will be generated for them.
            if frg_serie[0] not in ['y','b'] : continue


            for aa in sequence_heavy :
                if aa in labeling:
                    precursorMZ_heavy += labeling[aa].deltamass / z_parent

            frg_seq = sequence_heavy
            #b series
            if frg_serie == 'b':
                frg_seq = sequence_heavy[:frg_number]
            #y series
            if frg_serie == 'y':
                frg_seq = sequence_heavy[-frg_number:]

            for aa in frg_seq:
                if aa in labeling:
                    fragment_mz_heavy += labeling[aa].deltamass / frg_z

            #Check for C- and N-terminal labelings
            if 'C-term' in labeling:
                precursorMZ_heavy += labeling['C-term'].deltamass / z_parent
                if frg_serie == 'y':
                    fragment_mz_heavy += labeling['C-term'].deltamass / frg_z
            if 'N-term' in labeling:
                precursorMZ_heavy += labeling['N-term'].deltamass / z_parent
                if frg_serie == 'b':
                    fragment_mz_heavy += labeling['N-term'].deltamass / frg_z


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
            if lock : lock.acquire()
            writer.writerow(heavy_transition)
            time.sleep(0.001)
            if lock : lock.release()




if __name__ == '__main__':
    main(sys.argv[1:])
