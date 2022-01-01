#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Conversion plugins for various structure formats

The recognized structure formats are defined by subclassing StructureParser,
by convention these classes are named P_<format>.  The parser classes should
to override the parseLines() and toLines() methods of StructureParser.
Any structure parser needs to be registered in parser_index module.

For normal usage it should be sufficient to use the routines provided
in this module.

Content:
    StructureParser -- base class for concrete Parsers
    parser_index    -- dictionary of known structure formats
    getParser       -- factory for Parser at given format
    inputFormats    -- list of available input formats
    outputFormats   -- list of available output formats
"""

from diffpy.Structure import StructureFormatError
from diffpy.Structure.Parsers.structureparser import StructureParser
from diffpy.Structure.Parsers.parser_index_mod import parser_index

def getParser(format, **kw):
    """Return Parser instance for a given structure format.

    kw   -- keyword arguments passed to the Parser init function.

    Raises StructureFormatError exception when format is not defined.
    """
    if format not in parser_index:
        emsg = "no parser for '%s' format" % format
        raise StructureFormatError(emsg)
    pmod = parser_index[format]['module']
    print('CHEECKK', pmod)
    #pm = None
    import_cmd = 'from diffpy.Structure.Parsers import %s as pm' % pmod
    exec(import_cmd)
    from diffpy.Structure.Parsers import P_cif as pm
    return pm.getParser(**kw)

def inputFormats():
    """Return list of implemented input structure formats"""
    input_formats = [ fmt for fmt, prop in parser_index.items()
            if prop['has_input'] ]
    input_formats.sort()
    return input_formats

def outputFormats():
    """return list of implemented output structure formats"""
    output_formats = [ fmt for fmt, prop in parser_index.items()
            if prop['has_output'] ]
    output_formats.sort()
    return output_formats

# End of file
