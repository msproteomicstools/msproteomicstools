#!/usr/bin/python
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
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

from __future__ import print_function
from ..util.utils import *
import sys
import re
from ..data_structures.peak import Peak

class Library:
    """This class contains one spectral library, whatever that means.
    It provides an read/write interface to the database.
    It provides an read/write interface to the SpectraST *.splib and *.pepidx files.
    One can easily add spectra or retrive the spectra"""

    def __init__(self, lkey=None):
        self.seqHandler_hash = {}
        self.verbose = False
        if lkey: self.set_library_key( lkey )

        #self.pepidx = { sequence/mod : binindex }
        self.pepidx = {}

    def init_with_self(self, library):
        """Initialize with another library. Doesnt do a very deep copy """
        self.seqHandler_hash = library.seqHandler_hash
        self.set_library_key(library.library_key)
        for seq, handler in self.seqHandler_hash.iteritems():
            h = SequenceHandler()
            h.spectras = handler.spectras
            h.meta = handler.metas
            self.seqHandler_hash[ seq ] = h

    def __iter__(self):
        for handler in self.seqHandler_hash.values():
            for spectra in handler.spectras:
                yield spectra

    def nr_unique_peptides(self):
        return len( self.seqHandler_hash )

    def set_library_key(self, lkey ):
        self.library_key = lkey

    def measure_nr_spectra(self):
        leng = 0
        for key, h in self.seqHandler_hash.iteritems():
            leng += len( h.spectras )
        return leng

    def remove_duplicate_entries(self):
        #remove duplicates
        for handler in self.seqHandler_hash.values():
            handler.remove_duplicate_entries()
            if handler.empty():
                del self.seqHandler_hash[ s.sequence ]

    def get_spectra_by_sequence(self, sequence):
        """
        Get all spectra that match a specific sequence
        """
        result = []
        handler = self.seqHandler_hash[ sequence ]
        for spectra in handler.spectras:
            result.append(spectra)
        return result

    def get_all_spectra(self):
        result = []
        for handler in self.seqHandler_hash.values():
            for spectra in handler.spectras:
                result.append(spectra)
        return result

    def all_spectra(self):
        """
        Iterate over all specra in the library

        Yield:
            spectrum(:class:`.Spectra`): current spectrum
        """
        for handler in self.seqHandler_hash.values():
            for spectra in handler.spectras:
                yield spectra

    def annotate_with_libkey(self):
        """
        Annotate spectra with the key of the current library
        """
        for handler in self.seqHandler_hash.values():
            for s in handler.spectras:
                prot = s.comment.find( ' Protein=')
                length = s.comment[prot+1:].find( ' ')
                protein_str = s.comment[prot:prot+length]
                protein_desc = protein_str.split("=")[1]
                all_proteins = protein_desc.split("/")

                annotated = [ str(self.library_key) + "_" + p
                             for p in all_proteins[1:] ]
                annotated = [ str( len( annotated ) ) ] +  annotated
                annotated = reduce( lambda x,y: x + '/' + y, annotated )

                assert len(all_proteins) -1 == int( all_proteins[0] )
                s.comment = s.comment[:prot] + ' Protein=' + \
                        annotated + \
                        s.comment[prot+length:]

    def delete_reverse_spectra(self):
        for handler in self.seqHandler_hash.values():
            to_remove = []
            for s in handler.spectras:
                prot = s.comment.find( ' Protein=')
                length = s.comment[prot+1:].find( ' ')
                protein_str = s.comment[prot:prot+length]
                protein_desc = protein_str.split("=")[1]
                all_proteins = protein_desc.split("/")
                assert len(all_proteins) -1 == int( all_proteins[0] )
                rev = [r for r in [ p.find('reverse') for p in all_proteins]
                       if r != -1]
                if len(rev) == int( all_proteins[0]):
                    #only if all of them are reverse
                    to_remove.append( s)
            for s in to_remove:
                handler.remove( s )
            if handler.empty():
                del self.seqHandler_hash[ s.sequence ]

    def add_spectra(self, s):
        if self.seqHandler_hash.has_key(s.sequence):
            self.seqHandler_hash[ s.sequence ].add_spectra(s)
        else:
            self.seqHandler_hash[ s.sequence ] = SequenceHandler()
            self.seqHandler_hash[ s.sequence ].add_spectra(s)

    ##########################################
    #DB Functions

    def write_toDB(self, db, cursor):
        """
        Write all spectra into a SQL database
        """
        for handler in self.seqHandler_hash.values():
            for spectra in handler.spectras:
                spectra.save( db, cursor)

    def delete_library_from_DB( self, library_key, db):
        """
        Delete current library from SQL database
        """
        c = db.cursor()
        c.execute( "delete from hroest.specLibSpectra\
                  where id in (select id from \
                  hroest.specLibMeta where library_key = %s)" % library_key )
        query = "delete from hroest.specLibMeta\
                where library_key = %s" % library_key
        c.execute( query )

    def find_by_sequence(self, sequence, db):
        """
            This function can be used to access spectra using a sequence search
        """
        query = """
        SELECT id FROM hroest.specLibMeta
        WHERE sequence = '%s'""" % sequence
        self.find_by_sql( query, db)

    def find_by_sql_fast(self, subQuery, db, tmp_db):
        """
            This function can be used to access spectra using an sql query.
            The query should produce a single coloumn with spectra_keys (ids)
            which MUST be called tmp_spectra_keys.
            You need create table privileges in the databse tmp_db for this.
            But it can be 400x times faster than plain find_by_sql.
        """

        import os
        s = Spectra()
        pid = os.getpid()
        c = db.cursor()
        c.execute( "DROP TABLE IF EXISTS %s.%s_tmp_delete" % (tmp_db, pid))
        createTableQuery = "CREATE TABLE %s.%s_tmp_delete AS (%s)" % (tmp_db,
                pid, subQuery)
        #createTableQuery = "CREATE TEMPORARY TABLE %s.%s_tmp_delete\
        #       AS (%s)" % (tmp_db, pid, subQuery)
        c.execute(createTableQuery)
        c.execute( "ALTER TABLE %s.%s_tmp_delete ADD INDEX(tmp_spectra_keys)" %(
            tmp_db, pid))
        subQuery2 = "SELECT * FROM %s.%s_tmp_delete" % (tmp_db, pid)

        query = """
        SELECT %s FROM hroest.specLibMeta
        INNER JOIN hroest.specLibSpectra ON specLibMeta.id = specLibSpectra.id
        WHERE specLibMeta.id IN (%s)""" % (s.get_spectra_headers() + ", "+
                                           s.get_meta_headers(), subQuery2)
        query = query.replace( "id, ", "specLibMeta.id as id, " )

        t = utils.db_table(c)
        t.read( query )
        self._read_from_db_table( t)
        c.execute( "DROP TABLE %s.%s_tmp_delete" % (tmp_db, pid))

    def find_by_sql(self, query_in, db):
        """
            This function can be used to access spectra using an sql query.
            The query should produce a single coloumn with spectra_keys.
            This can be very slow, use find_by_sql_fast instead (~400x faster).
        """

        query = """
        SELECT * FROM hroest.specLibMeta #USE INDEX PRIMARY
        INNER JOIN hroest.specLibSpectra ON specLibMeta.id = specLibSpectra.id
        WHERE specLibMeta.id IN (%s)""" % query_in
        c = db.cursor()
        t = utils.db_table(c)
        import time
        start = time.time()
        t.read( query )
        end = time.time()
        print("query took %s s" % ( end - start))
        self._read_from_db_table( t)

    def read_fromDB(self, library_key, db):
        """This function can be used to access one complete library from the DB."""
        s = Spectra()
        query = """
        SELECT %s FROM hroest.specLibMeta
        INNER JOIN hroest.specLibSpectra ON specLibMeta.id = specLibSpectra.id
        WHERE library_key = %s""" % (s.get_spectra_headers() + ", "+
                                     s.get_meta_headers(), library_key)
        query = query.replace( "id, ", "specLibMeta.id as id, " )
        c = db.cursor()
        t = utils.db_table(c)
        t.read( query )
        self._read_from_db_table(t)

    @staticmethod
    def read_from_db_to_file(library_key, db, filePrefix):
        """This function can be used to access one complete library from the
        DB directly to a file."""
        s = Spectra()
        query = """
        SELECT %s FROM hroest.specLibMeta
        INNER JOIN hroest.specLibSpectra ON specLibMeta.id = specLibSpectra.id
        WHERE library_key = %s""" % (s.get_spectra_headers() + ", "+
                                     s.get_meta_headers(), library_key)
        query = query.replace( "id, ", "specLibMeta.id as id, " )
        c = db.cursor()
        t = utils.db_table(c)
        t.read( query )
        generator = Library._get_spectra_generator_from_dbTable(t)
        Library._write_to_file( generator , filePrefix)

    def _read_from_db_table(self, t):
        for row in t.rows():
            s = Spectra()
            s._initialize_dbtable(row, t, s.meta_cols )
            s._initialize_dbtable(row, t, s.spectra_cols )
            self.add_spectra( s )

    @staticmethod
    def _get_spectra_generator_from_dbTable(t):

        def myIterator():
            for row in t.rows():
                s = Spectra()
                s._initialize_dbtable(row, t, s.meta_cols )
                s._initialize_dbtable(row, t, s.spectra_cols )
                yield s

        return myIterator()

    ##########################################
    #File Functions

    def write(self, filePrefix, append=False):
        """Write the current library to a file."""

        Library._write_to_file( self , filePrefix, append)

    @staticmethod
    def _write_to_file(spectra_iterator, filePrefix, append):
        """Write a library to a file."""
        splibFile = filePrefix + '.splib'
        pepIdxFile = filePrefix + '.pepidx'
        if append: mode = 'a'
        else: mode = 'w'
        fs = open( splibFile,   mode )
        fp = open( pepIdxFile,  mode )
        i = 0
        for spectra in spectra_iterator:
            spectra.binindex = fs.tell()
            spectra.LibID = i
            #very important to add a final new line
            fs.write( spectra.to_splib_str() + '\n')
            fp.write( spectra.to_pepidx_str() )
            i += 1

        fs.close()
        fp.close()

    def write_sorted(self, filePrefix):
        splibFile = filePrefix + '.splib'
        pepIdxFile = filePrefix + '.pepidx'
        fs = open( splibFile, 'w')
        fp = open( pepIdxFile, 'w')
        i = 0
        for sequence in sorted(self.seqHandler_hash):
            handler = self.seqHandler_hash[ sequence ]
            for spectra in handler:
                spectra.binindex = fs.tell()
                spectra.LibID = i
                #very important to add a final new line
                fs.write( spectra.to_splib_str() + '\n')
                fp.write( spectra.to_pepidx_str() )
                i += 1
        fs.close()
        fp.close()

    def get_first_offset(self, splibFileName) :

        line_offset = 0
        last_offset = 0

        fs = open (splibFileName , 'r')

        while True:
            line_offset = fs.tell()
            row = fs.readline()

            #~ print row
            if len(row) == 0 : break
            if len(row) > 5 and row[:5] == 'Name:' : break

        last_offset = line_offset

        #cleanup
        fs.close()

        return last_offset

    def get_fileheader(self, splibFileName):
        """Get the header preceding the first spectrum in a spectrast file."""
        
        fs = open( splibFileName , 'r' )
        
        firstRow = fs.readline()

        stack = [ firstRow ]
        while True:
            line_offset = fs.tell()
            row = fs.readline()

            #~ print row
            if len(row) == 0:
                line_offset = offset
                break
            if len(row) > 5 and row[:5] == 'Name:': break
            stack.append( row )

        #cleanup
        fs.close()
        
        return stack
    
    def get_rawspectrum_with_offset(self, splibFileName, offset) :
        """Get a raw spectrum as it is from a spectrast file by using an offset to locate it."""
        
        assert splibFileName[-6:-3] == '.sp'
        
        fs = open( splibFileName , 'r' )
        
        binindex = int ( offset )
        
        fs.seek( binindex )
        firstRow = fs.readline()

        stack = [ firstRow ]
        while True:
            line_offset = fs.tell()
            row = fs.readline()

            #~ print row
            if len(row) == 0:
                line_offset = offset
                break
            if len(row) > 5 and row[:5] == 'Name:': break
            stack.append( row )

        #cleanup
        fs.close()
        
        return stack
        
    def read_sptxt_with_offset(self, splibFileName, offset) :
        """Read a sptxt spectra library file by using an offset to keep memory free"""

        assert splibFileName[-6:-3] == '.sp'

        line_offset = 0
        last_offset = 0

        fs = open( splibFileName , 'r' )

        binindex = int ( offset )

        s = Spectra()
        ####################################
        #go to binary offset in the splib file, read spectrum
        fs.seek( binindex )
        firstRow = fs.readline()
        #~ assert firstRow[:5] == 'Name:', error_notmatch

        stack = [ firstRow ]
        while True:
            line_offset = fs.tell()
            row = fs.readline()

            #~ print row
            if len(row) == 0:
                line_offset = offset
                break
            if len(row) > 5 and row[:5] == 'Name:': break
            stack.append( row )

        s.parse_sptxt( stack )
        last_offset = line_offset

        #cleanup
        fs.close()

        return last_offset, s

    def read_spectrum_sptxt_idx(self, splibFileName, idx ,library_key) :
        """"Fetch a spectrum from the spectral library, by using the binary index"""

        assert splibFileName[-6:] == '.sptxt'

        fs = open( splibFileName , 'r' )
        error_notmatch = \
    """\n    An error occured while reading: \n    %s.
    The pepidx file does not correspond to splib file. Did they come from a
    different library or is the splib file in binary format? If so, try to
    run SpectraST in library creation mode with the option -c_BIN! set.
        """ % splibFileName

        binindex = int(idx)

        s = Spectra()
        ####################################
        #go to binary offset in the splib file, read spectrum
        fs.seek( binindex )
        firstRow = fs.readline()
        assert firstRow[:5] == 'Name:', error_notmatch
        stack = [ firstRow ]
        while True:
            row = fs.readline()
            print(row)
            if len(row) == 0: break
            if len(row) > 5 and row[:5] == 'Name:': break
            stack.append( row )
        s.parse_sptxt( stack )

        #cleanup
        fs.close()

        return s


    def read_sptxt_pepidx(self, splibFileName, pepidxFileName,library_key):
        """Read directly from a spectral library into memory. """

        assert splibFileName[-6:] == '.splib'
        assert pepidxFileName[-7:] == '.pepidx'
        fs = open( splibFileName , 'r')
        fp = open( pepidxFileName , 'r')
        error_notmatch = \
    """\n    An error occured while reading: \n    %s, %s
    The pepidx file does not correspond to splib file. Did they come from a
    different library or is the splib file in binary format? If so, try to
    run SpectraST in library creation mode with the option -c_BIN! set.
        """ % (splibFileName, pepidxFileName)
        while True:
            row = fp.readline()
            if len(row) == 0: break
            if row[0] == '#': continue
            srow = row.split("\t")
            sequence = srow[0]
            modifications = srow[1]

            #there may be two or more binary indices per line
            binindex_all = srow[2]
            for b in binindex_all.split():
                binindex = int(b)
                s = Spectra()
                ####################################
                #go to binary offset in the splib file, read spectrum
                fs.seek( binindex )
                firstRow = fs.readline()
                assert firstRow[:5] == 'Name:', error_notmatch
                stack = [ firstRow ]
                while True:
                    row = fs.readline()
                    if len(row) == 0: break
                    if len(row) > 5 and row[:5] == 'Name:': break
                    stack.append( row )
                s.parse_sptxt( stack )
                s.add_meta( sequence, modifications, library_key)
                #s.set_library_key( library_key )
                self.add_spectra( s )

        #cleanup
        fs.close()
        fp.close()


    def read_sptxt(self, filename):
        #raise ValueError, "deprecated"
        f = open( filename, 'r')
        line_offset = 0
        last_offset = 0
        i = 0
        self.seqHandler_hash = {}
        stack = []
        while True:
            line_offset = f.tell()
            row = f.readline()
            if len(row) > 0 and row[0] == '#': continue
            if len(row) == 0 or (row[:5] == 'Name:' and len(stack) > 0):
                i += 1
                s = Spectra()
                s.parse_sptxt( stack )
                #these do not necessarely agree with the values in *.pepidx
                #since those COULD be built from the binary file.
                s.binindex = last_offset
                self.add_spectra( s )
                stack=[]
                last_offset = line_offset
                if len(row) == 0: break
            stack.append( row )
        f.close()

    def read_pepidx(self, filename):
        #~ raise ValueError, "deprecated"
        f = open( filename , 'r')
        pepidx = f.readlines()
        f.close()
        for row in pepidx:
            if row[0] == '#': continue
            srow = row.split("\t")
            sequence = srow[0]
            modifications = srow[1]
            binindex = srow[2]

            #self.pepidx = { sequence/mod : binindex }
            sequencemod = sequence + '/' + modifications
            self.pepidx[sequencemod] = binindex

            #~ meta = SpectraMeta()
            #~ meta.initialize(sequence, modifications, binindex, self.library_key)
            #~ #TODO catch error ==> hopefully no error,
            #~ #we should find all the spectra again
            #~ handler = self.seqHandler_hash[ sequence ]
            #~ handler.add_meta( meta )

    @staticmethod
    def read_library_to_db(splibFileName, pepidxFileName, db, library_key):
        """Read directly from a spectral library into the database. """
        assert splibFileName[-6:] == '.splib'
        assert pepidxFileName[-7:] == '.pepidx'
        fs = open( splibFileName , 'r')
        fp = open( pepidxFileName , 'r')
        while True:
            row = fp.readline()
            if len(row) == 0: break
            if row[0] == '#': continue
            srow = row.split("\t")
            sequence = srow[0]
            modifications = srow[1]

            #there may be two or more binary indices per line
            binindex_all = srow[2]
            for b in binindex_all.split():
                binindex = int(b)
                s = Spectra()
                ####################################
                #go to binary offset in the splib file, read spectrum and save
                fs.seek( binindex )
                firstRow = fs.readline()
                assert firstRow[:5] == 'Name:'
                stack = [ firstRow ]
                while True:
                    row = fs.readline()
                    if len(row) == 0: break
                    if len(row) > 5 and row[:5] == 'Name:': break
                    stack.append( row )
                s.parse_sptxt( stack )
                s.add_meta( sequence, modifications, library_key)
                s.save( db )

        #cleanup
        fs.close()
        fp.close()

    def count_modifications(self):
        myhash = {}
        for handler in self.seqHandler_hash.values():
            #for meta in handler.metas:
            for meta in handler.spectras:
                spl = meta.ptm_string.split("|")
                mod = spl[1].split( "/" )
                nr_mod = int(mod[0])
                assert nr_mod == len( mod ) -1
                mods = mod[1:]
                for m in mods:
                    values = m.split(",")
                    #['3', 'C', 'Propionamide']
                    myhash[ values[2] ] = 0

class SequenceHandler:
    """
    Container class of spectra with the same sequence in a spectral library

    Acts as a container of all spectra mapping to the same sequence inside a spectral library
    """

    def __init__(self):
        self.spectras = []
        self.metas = []

    def remove_duplicate_entries(self):
        ptm_dic = []
        res = []
        for spectra in self.spectras:
            if spectra.ptm_string not in ptm_dic:
                res.append( spectra )
            ptm_dic.append( spectra.ptm_string)
        self.spectras = res

    def init_with_self(self, handler):
        for spectra in handler.spectras:
            s = Spectra()
            s.init_with_self( spectra )
            self.spectras.append( s )
        for meta in handler.metas:
            meta = SpectraMeta()
            meta.init_with_self( meta )
            self.metas.append( meta )

    def add_spectra(self, spectra):
        self.spectras.append( spectra )

    def add_spectra_no_duplicates(self, spectra):
        for s in self.spectras:
            if (s.meta.modifications == spectra.meta.modifications and
                s.LibID == spectra.LibID and s.full == spectra.full): return
        self.spectras.append( spectra )

    def add_meta( self, meta):
        self.metas.append( meta)

    def remove(self, s):
        assert s in self.spectras
        self.spectras.remove( s)
        assert s not in self.spectras

    def empty(self):
        return len(self.spectras) == 0

class Spectra:
    """
    A single spectrum inside a spectral library
    
    """

    def __init__(self):
        self.initialize()

    def initialize(self):
        """
        Initialize spectrum
        """
        self.compress_spectra = ''
        self.meta = None
        self.phosphos = []
        self.methyls = []
        self.acetyls = []
        self.carbas = []
        self.oxidations = []
        self.icat = []
        self.other_known = []
        self.other = []
        self.meta_cols = {
            'library_key'           :   -999,
            'sequence'              :   1,
            'charge'                :   0,
            'number_peaks'          :   0,
            'MW'                    :   0,
            'ConsFracAssignedPeaks' :   0,
            'OrigMaxIntensity'      :   0,
            'precursorMZ'           :   0,
            'nr_replicates'         :   0,
            'phospho_len'           :   2,
            'methyl_len'            :   2,
            'acetyl_len'            :   2,
            'carbamido_len'         :   2,
            'oxidations_len'        :   2,
            'icat_len'              :   2,
            'other_known_len'       :   2,
            'other_len'             :   2
        }
        #optional
        self.peptide = ''
        self.protein_ac = ''
        self.OrigMaxIntensity = -1
        self.ConsFracAssignedPeaks = -1
        self.RetTime = -1
        self.RetTime_detected = False
        self.iRT = -1
        self.iRT_detected = False
        self.spectra_cols = {
            'id'                :   0,
            'status'            :   1,
            'name'              :   1,
            'full_name'         :   1,
            'comment'           :   1,
            'ptm_string'        :   1,
            'peptide'           :   1,
            'protein_ac'        :   1,
            'compress_spectra'  :   5
        }

    #TODO use @property for the following
    #
    def phospho_len(self):
        return len(self.phosphos)
    #
    def methyl_len(self):
        return len(self.methyls)
    #
    def acetyl_len(self):
        return len(self.acetyls)
    #
    def carbamido_len(self):
        return len(self.carbas)
    #
    def oxidations_len(self):
        return len(self.oxidations)
    #
    def icat_len(self):
        return len(self.icat)
    #
    def other_known_len(self):
        return len(self.other_known)
    #
    def other_len(self):
        return len(self.other)
    #
    def _initialize_dbtable(self, row, t, hash):
        for attr, value in hash.iteritems():
            exec("self.%s = t.row( row, '%s')" % (attr, attr))

    def to_pepidx_str(self):
        """
        Convert spectrum object to pepidx format
        """
        s = ''
        s += self.sequence + '\t'
        s += self.ptm_string + '\t'
        s += repr(self.binindex) + '\n'
        return s

    def to_splib_str(self):
        """
        Convert spectrum object to splib format
        """
        s = ''
        s += 'Name: '+              self.name + '\n'
        s += 'LibID: '+        str( self.LibID ) + '\n'
        s += 'MW: '+           str( self.MW ) + '\n'
        s += 'PrecursorMZ: '+  str( self.precursorMZ) + '\n'
        s += 'Status: '+            self.status + '\n'
        s += 'FullName: '+          self.full_name + '\n'
        s += 'Comment: '+           self.comment + '\n'
        s += 'NumPeaks: '+     str( self.number_peaks) + '\n'
        s += self.compress_spectra
        return s

    #def init_with_self(self, spectra):
    #    import copy
    #    self.name =        copy.deepcopy(spectra.name)
    #    self.libId =       copy.deepcopy(spectra.libId)
    #    self.mW =          copy.deepcopy(spectra.mW)
    #    self.precursor =   copy.deepcopy(spectra.precursor)
    #    self.status =      copy.deepcopy(spectra.status)
    #    self.full =        copy.deepcopy(spectra.full)
    #    self.comment =     copy.deepcopy(spectra.comment)
    #    self.nPeaks =      copy.deepcopy(spectra.nPeaks)
    #    self.peptide =     copy.deepcopy(spectra.peptide)
    #    self.protein_as =  copy.deepcopy(spectra.protein_as)
    #    for peak in spectra.peaks:
    #        p = Peak()
    #        p.init_with_self( peak )
    #        self.peaks.append( p )

    def get_known_modifications(self):
        known = [
        'Phopsho'
        'Acetyl',
        'Methyl',
        'Carbamidomethyl',
        'Oxidation',
        'Gln->pyro-Glu',
        'ICAT_light',
        'ICAT_heavy',
        'Pyro-cmC',
        'AB_old_ICATd8',
        'AB_old_ICATd0'
        ]
        return known

    def analyse_mod(self):
        #the ptm_string looks like this
        #
        # | 2|2/7,S,Phospho/-2,K,Methyl| |
        #
        # split by "|" and you get: charge, modifications
        # split modifications by "/" and you get:
        #   -- number of modifications
        #   -- modifications themselves
        #
        # then split each mod_string by "," and you get: position, aa, type
        #
        spl =  self.modifications.split("|")
        self.charge = int( spl[0])
        mod = spl[1].split( "/" )
        nr_mod = int(mod[0])
        assert nr_mod == len( mod ) -1
        mods = mod[1:]
        known = self.get_known_modifications()
        for m in mods:
            values = m.split(",")
            if values[2] == 'Phospho': self.phosphos.append([values[0], values[1]] )
            elif values[2] == 'Methyl': self.methyls.append([values[0], values[1]] )
            elif values[2] == 'Acetyl': self.acetyls.append([values[0], values[1]] )
            elif values[2] == 'Carbamidomethyl': self.carbas.append([values[0], values[1]] )
            elif values[2] == 'Oxidation': self.oxidations.append([values[0], values[1]] )
            elif values[2] in ['ICAT_heavy', 'ICAT_light']: self.icat.append([values[0], values[1]] )
            elif values[2] in known: self.other_known.append([ values[0], values[1]])
            else: self.other.append([ values[0], values[1] ])

    def validate(self):
        spl =  self.modifications.split("|")
        mod = spl[1].split( "/" )
        nr_mod = int(mod[0])
        assert nr_mod == (
            len(self.phosphos) +
            len(self.methyls) +
            len(self.acetyls) +
            len(self.carbas) +
            len(self.oxidations) +
            len(self.icat) +
            len(self.other_known) +
            len(self.other)
        )

        #import DDB
        #import silver
        #R = silver.Residues.Residues('mono')
        #peptide = DDB.Peptide()
        #peptide.set_sequence( self.meta.sequence )
        #S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        #S.ion_charge = s.charge
        #S.construct_from_peptide( peptide.get_modified_sequence('SEQUEST'),
        #                         R.residues, R.res_pairs)
        #TODO turn on this assertion
        #assert S.peptide_mass / S.ion_charge - self.MW < 1e-2

    def parse_sptxt(self, stack):
        """
        Parse an sptxt entry and initialize spectrum
        """

        def median(y):
            z = len(y)
            if not z%2: return (y[(z/2)-1] + y[z/2]) / 2
            else: return y[z/2]


        self.initialize()
        peaks_start = False
        peaks = []
        for row in stack:
            if row[:5] == 'Name:':            self.name = row[6:-1]
            elif row[:6] == 'LibID:':         self.LibID = int(row[7:-1])
            elif row[:3] == 'MW:':            self.MW = float(row[4:-1])
            elif row[:12] == 'PrecursorMZ:':  self.precursorMZ = float(row[13:-1])
            elif row[:7] == 'Status:':        self.status = row[8:-1]
            elif row[:9] == 'FullName:':      self.full_name = row[10:-1]
            elif row[:8] == 'Comment:':       self.comment = row[9:-1]
            elif row[:9] == 'NumPeaks:':      self.number_peaks= int(row[10:-1]); peaks_start = True
            elif peaks_start and not (row.strip()  == ''):
                self.compress_spectra += row

        import re
        #TODO make only the 20 aa possible
        self.sequence =  re.sub( '[^A-Z]', '', self.name )
        mm = re.search( 'Nreps=(.*?)\s', self.comment )
        if mm: self.nr_replicates = int( mm.group(1).split("/")[0] )
        m = re.search( 'Pep=(.*?)\s', self.comment )
        if m: self.peptide = m.group(1)
        mm = re.search( 'Protein=(.*?)\s', self.comment )
        if mm: self.protein_ac = mm.group(1)
        mm = re.search( 'ConsFracAssignedPeaks=(.*?)\s', self.comment )
        if mm: self.ConsFracAssignedPeaks = float( mm.group(1) )
        mm = re.search( 'OrigMaxIntensity=(.*?)\s', self.comment )
        if mm: self.OrigMaxIntensity = float( mm.group(1) )
        mm = re.search( 'RetentionTime=(.*?)\s', self.comment )
        if mm:
            if ',' in mm.group(1)[:] :
                mm_split = mm.group(1).split(',')
                mm_split_float = [float(f) for f in mm_split]
                self.RetTime = mm_split_float[0]
                self.RetTime_detected = True
        mm = re.search( 'iRT=(.*?)$', self.comment )
        if mm:
            if ',' in mm.group(1)[:] :
                mm_split = mm.group(1).split(',')
                mm_split_float = [float(f) for f in mm_split]
                self.iRT = mm_split_float[0]
                self.iRT_detected = True

        try :
            self.comment_parsed = self.parse_comments(self.comment)
        except : pass
        if 'Se' in self.comment_parsed :
            self.searchEngineInfo = self.parse_SearchEngineInfo(self.comment_parsed['Se'])

    def add_meta(self, sequence, modifications, library_key):
        if self.sequence: assert sequence == self.sequence
        self.modifications = modifications
        self.library_key = library_key
        self.analyse_mod()

    def phosphos_len(self):
        #TODO dont use re
        import re
        ser = re.findall( '167', self.name)
        thr = re.findall( '181', self.name)
        tyr = re.findall( '243', self.name)
        return len(ser) + len(thr) + len(tyr)

    def _headers(self, hash):
        query = ''
        for attr, value in hash.iteritems():
            if value == 5:
                query += 'uncompress( %s ) as %s, ' % (attr, attr)
            else:
                query += attr + ', '
        return query[:-2]

    def get_meta_headers(self):
        return self._headers( self.meta_cols )

    def get_spectra_headers(self):
        return self._headers( self.spectra_cols )

    def find(self, id, db):
        query = 'SELECT ' + self.get_meta_headers()
        query += " FROM hroest.specLibMeta WHERE id = %s" % id
        c = db.cursor()
        t = utils.db_table(c)
        t.read( query )
        #there should be only one result
        assert t.result == 1
        row = t.fetchone()
        self._initialize_dbtable(row, t, self.meta_cols)
        #
        query = 'SELECT ' + self.get_spectra_headers()
        query += " FROM hroest.specLibSpectra WHERE id = %s" % id
        t = utils.db_table(c)
        t.read( query )
        #there should be only one result
        assert t.result == 1
        row = t.fetchone()
        self._initialize_dbtable(row, t, self.spectra_cols)

    def _meta_save(self, db, cursor):
        head = """INSERT INTO hroest.specLibMeta ("""
        query = '('
        for attr, value in self.meta_cols.iteritems():
            head += attr + ", "
            if value == 0:
                exec( "query +=str(self.%s) + ', '" % attr)
            if value == 1:
                exec( "query +='\\'%%s\\'' %% self.%s + ', '" % attr)
            if value == 2:
                exec( "query +=str(self.%s()) + ', ' " % attr)
        head = head[:-2] + ") VALUES "
        query = head + query[:-2] + ')'
        self.validate()
        cursor.execute(query)
        last_id = db.insert_id()
        return last_id

    def save(self, db):
        #first we get the next id back from the meta table
        #since id is contained in spectra_cols, we will also store it
        cursor = db.cursor()
        self.id = self._meta_save( db, cursor )
        self.comment = self.escape_string( self.comment )
        self.ptm_string = self.modifications
        head = """INSERT INTO hroest.specLibSpectra ("""
        query = '('
        for attr, value in self.spectra_cols.iteritems():
            head += attr + ", "
            if value == 0:
                exec( "query +=str(self.%s) + ', '" % attr)
            if value == 1:
                exec( "query +='\\'%%s\\'' %% self.%s + ', '" % attr)
            if value == 2:
                exec( "query +=str(self.%s()) + ', ' " % attr)
            if value == 5:
                exec( "query +='compress(\\'%%s\\')' %% self.%s + ', '" % attr)
        head = head[:-2] + ") VALUES "
        query = head + query[:-2] + ')'
        self.validate()
        cursor.execute(query)

    def escape_string(self, string):
        #TODO dont use re
        import re
        res = re.sub( '\'', "\\\'", string)
        return res

    def get_peaks(self):
        return [Peak(p, spectraST = True) for p in self.compress_spectra.split('\n') if not p == '']

    def is_tryptic(self):
        return self.comment.find( ' Pep=Tryptic ' ) != -1

    def parse_comments(self, comment):
        comment2 = ' ' + comment + ' dummy='
        fields = re.split(' (\w*)=(.*?) (\w*)=', comment2)
        fields = fields[1:-2]
        comment_fields = {}
        for i in range( 0, len(fields), 2):
           comment_fields[ fields[i] ] = fields[i+1]
        return comment_fields

    def parse_SearchEngineInfo(self, searchEngineInfo) :
        # 'Se': '2^M20:ex=4.5000e-04/0.0000e+00,fv=4.3803/1.3592,hs=29.5570/12.3713,is=26.0340/0.1320,pb=1.0000/0.0000,sc=75.3035/13.5632,sr=45.7465/21.3761^S19:dc=0.3386/0.0656,fv=3.7786/1.1032,pb=0.9970/0.0106,xc=2.7534/0.4410'

        #search engines
        engines = searchEngineInfo.split('^')[1:]
        numofengines = len(engines)

        searchEngineInfo_parsed = {}
        for engine in engines :
            eng_info = engine.split(':')
            if len(eng_info) > 1 :
                engine = eng_info[0][0]
                info   = eng_info[1].split(',')
            else :
                engine = 'Default'
                info   = eng_info.split(',')

            info_fields = {}
            for inf in info :
                field = inf.split('=')[0]
                value = inf.split('=')[1]
                sd = 'NA'
                if '/' in value :
                    sd = value.split('/')[1]
                    value = value.split('/')[0]
                info_fields[field] = (value,sd)
            searchEngineInfo_parsed[engine] = info_fields

        return searchEngineInfo_parsed

class Peak_old:
    """DEPRECATED: Represents one peak of a spectrum."""
    def __init__(self, str=None):
        if str: self.parse_str( str )

    def initialize(self, peak, intensity, peak_annotation, statistics):
        self.peak            = float(peak)
        self.intensity       = float(intensity)
        self.peak_annotation = peak_annotation
        self.statistics      = statistics
        self._parse_statistics( self.statistics )

    def init_with_self( self, peak):
        import copy
        self.peak            = copy.deepcopy(peak.peak)
        self.intensity       = copy.deepcopy(peak.intensity)
        self.peak_annotation = copy.deepcopy(peak.peak_annotation)
        self.statistics      = copy.deepcopy(peak.statistics)
        self._parse_statistics( self.statistics )

    def _parse_statistics(self, stat):
        if stat.strip() == '':
            self.nr_replicates = self.nr_replicates_necessary = -1
            self.mzStdDev = self.intensityStdDevOverMean = -1
            return
        stripped = stat.strip()
        split = stripped.split(" ")
        s_lib_occurence = split[0].split("/")
        self.nr_replicates = int(s_lib_occurence[0])
        self.nr_replicates_necessary = int(s_lib_occurence[1])
        statistics = split[1].split("|")
        if len(statistics) > 0: self.mzStdDev = float(statistics[0])
        else: self.mzStdDev = -1
        if len(statistics) > 1: self.intensityStdDevOverMean = float(statistics[1])
        else: self.intensityStdDevOverMean = -1

    def parse_str(self, peak):
        speak = peak.split( "\t")
        peak = float( speak[0] )
        intensity = float( speak[1] )
        if len(speak) > 2: peak_annotation = speak[2]
        if len(speak) > 3: statistics = speak[3]
        self.initialize(peak, intensity, peak_annotation, statistics)

    def to_write_string(self):
        s = ''
        s += repr(self.peak) + '\t'
        s += repr(self.intensity) + '\t'
        s += (self.peak_annotation) + '\t'
        s += (self.statistics) + '\t'
        return s

    def save(self, db, cursor, parent_key):
        head = """INSERT INTO hroest.specLibPeak (
            spectra_key, mz, intensity, peak_annotation, statistics,
            nr_replicates, mzStdDev, intensityStdDevOverMean)
            VALUES \n"""
        query = head +"""(%s, %s, %s, '%s', '%s',
                          %s, %s, %s) """ %\
        (parent_key, self.peak, self.intensity, self.peak_annotation,
                                                                self.statistics,
         self.nr_replicates, self.mzStdDev, self.intensityStdDevOverMean)

        cursor.execute(query)

