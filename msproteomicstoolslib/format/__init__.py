__all__ = ["pepXMLReader","methodDamReader","methodMethReader","mzXMLreader", "ProteinDB","speclib_db_lib"]
import sys
import csv

maxInt = sys.maxsize
decrement = True
while decrement:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.
    # http://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True
