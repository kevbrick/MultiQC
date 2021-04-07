
import numpy as np
from collections import OrderedDict

_SEP = '\t'
_KEY_SEP = '/'

# some good-looking presets for pair-type report section:
allelictypes_common = OrderedDict([('genome1'         , '#33a02c'), # Green
                                   ('genome2'         , '#6a3d9a'), # Violet
                                   ('unassignedN'     , '#ff7f00'), # orange
                                   ('unassigned_other', '#fff900'), # Yellow
                                   ('conflicting'     , '#e31a1c'), # Grey
                                   ('other'           , '#cccccc')])

class ParseError(Exception):
    pass

def read_allelicstatus_stats(file_handle):
    """create instance of PairCounter from file
    Parameters
    ----------
    file_handle: file handle
    Returns
    -------
    PairCounter
        new PairCounter filled with the contents of the input file
    """
    # fill in from file - file_handle:
    stat_from_file = OrderedDict()

    stat_from_file['allelic_status'] = OrderedDict([
        ('genome1', 0),
        ('genome2', 0),
        ('unassignedN', 0),
        ('unassigned_other', 0),
        ('conflicting', 0),
        ('other', 0)
        ])
    #
    # stat_from_file['genome1'] = 0
    # stat_from_file['genome2'] = 0
    # stat_from_file['unassigned_N'] = 0
    # stat_from_file['unassigned_other'] = 0
    # stat_from_file['conflicting'] = 0
    # stat_from_file['other'] = 0
    for l in file_handle:
        fields = l.strip().split(_SEP)
        if len(fields) == 0:
            # skip empty lines:
            continue
        if len(fields) != 2:
            # expect two _SEP separated values per line:
            raise ParseError(
                '{} is not a valid stats file'.format(file_handle.name))
        # extract key and value, then split the key:
        putative_key, putative_val =  fields[0], fields[1]
        key_fields = putative_key.split(_KEY_SEP)
        # we should impose a rigid structure of .stats or redo it:

        if len(key_fields) == 2:
            key = key_fields.pop(0)
            stat_from_file[key][key_fields[0]] = int(fields[1])
            print(key_fields[0] + " - " + fields[1] + "\n")
        else:
            raise ParseError(
                '{} is not a valid stats file: {} section implies 2 identifiers'.format(file_handle.name,key))

    #stat_from_file['cis_percent'] = stat_from_file['cis']/stat_from_file['total_nodups']*100.

    #
    return stat_from_file
