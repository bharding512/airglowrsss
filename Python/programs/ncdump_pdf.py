# This script reads a netCDF4 file and outputs a formatted PDF with a description of the 
# data product, dimensions, and variables. All information comes from metadata in the 
# netCDF4 file. Only Python 2.7 is supported. To run this script:
#
# $ python ncdump_pdf.py /path/to/ICON_L2_xxxx.NC
#
# A pdf will be created in the current working directory, with a name that is generated
# from the NC filename.

import numpy as np
import reportlab
import netCDF4
import datetime
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
import sys

PARAG_LIMIT = 17 # Limit the number of paragraphs in Var_Notes, otherwise reportlab will crash.

##### Define input and output filenames
usagestr = 'usage: python ncdump_pdf.py /path/to/ICON_L2_xxxx.NC'
if len(sys.argv)!=2:
    print usagestr
    sys.exit(1)
fn_in = sys.argv[1]

stub = '_'.join(fn_in.split('/')[-1].split('_')[:4])
vers = fn_in.split('/')[-1].split('_')[-1][:3]
fn_out = '%s_%s.pdf' % (stub, vers) # in current directory


# Helper function to determine if the attribute is filled in or not
def isvalid(x):
    '''
    Return True if the attribute is OK, meaning it is a:
    - string
    - multi-length string
    - something else that isn't NaN
    '''
    # isnan(s) results in an error if s is a string, so:
    if isinstance(x,(str, unicode)):
        return True
    if hasattr(x, '__len__'):
        return True
    else:
        return ~np.isnan(x) 

#################################################################################
######################## Extract information from file ##########################
#################################################################################

# Grab important global attributes, and all variables.
# Organize variables by Var_Type, and grab all important variable attributes
try:
    d = netCDF4.Dataset(fn_in)
except:
    raise IOError('File not found: %s' % fn_in)
    

a = {} # all global attributes
v = {} # keys are "data", "support_data", etc. Whatever the Var_Type can be (not case sensitive)
dims = {} # all dimensions
error_notes = [] # Notes to be printed in red at the top

try: # if something happens, make sure to close file
    
    for attr in d.ncattrs():
        a[attr.lower()] = d.getncattr(attr) # convert attribute names to lower case
        
    for dim in d.dimensions:
        dims[dim] = d.dimensions[dim].size
        if d.dimensions[dim].isunlimited():
            dims[dim] = 'unlimited'
    
    for var_name in d.variables:
        var = d[var_name]
        
        attrs = var.ncattrs()
        
        # Read variable attributes from netCDF file and convert to dictionary
        var_dict = {}
        var_dict['name'] = var_name
        for attr in var.ncattrs():
            var_dict[attr.lower()] = var.getncattr(attr) # convert attribute names to lower case
        # Read dimensions
        var_dict['dims'] = var.dimensions
        
        # Special case for error checking (most error checking is done below)
        if 'var_type' not in var_dict.keys() or not isvalid(var_dict['var_type']):
            var_dict['var_type'] = 'MISSING Var_Type'
            error_notes.append('ERROR: Missing variables attribute %s : Var_Type'%(var_dict['name']))

        var_type = var_dict['var_type'].lower() # not case sensitive, as per Tori's request
        if var_type not in ['ignore_data']:
            if var_type not in v.keys():
                v[var_type] = []

            v[var_type].append(var_dict)
    
    d.close()
except: 
    d.close()
    raise
    
#### Global attributes which are required
req_a = ['text_supplement','description','data_type','software_version']
for attr in req_a:
    if attr not in a.keys() or not isvalid(a[attr]):
        a[attr] = 'MISSING %s' % attr
        error_notes.append('ERROR: Missing global attribute "%s"'%attr)

#### Populate variable attributes which are required, and fill in those that are optional with blanks
req_v = ['var_type', 'var_notes', 'catdesc']
filler_v = ['units']
for var_type in v:
    for var_dict in v[var_type]:
        for attr in req_v:
            if attr not in var_dict.keys() or not isvalid(var_dict[attr]):
                var_dict[attr] = 'MISSING %s' % attr
                error_notes.append('ERROR: Missing variables attribute %s : %s'%(var_dict['name'],attr))
        for attr in filler_v:
            if attr not in var_dict.keys() or not isvalid(var_dict[attr]):
                # No need to throw error
                var_dict[attr] = ''
                
#### If Units are given as a list, convert it to a string
str_v = ['units']
for var_type in v:
    for var_dict in v[var_type]:
        for attr in str_v:
            if isinstance(var_dict[attr],list):
                var_dict[attr] = ', '.join(var_dict[attr])


#### Convert Var_Notes and Text_Supplement to multi-strings, if they are not already.
if isinstance(a['text_supplement'], (str, unicode)):
    a['text_supplement'] = [a['text_supplement']]
for var_type in v:
    for var_dict in v[var_type]:
        if isinstance(var_dict['var_notes'], (str,unicode)):
            var_dict['var_notes'] = [var_dict['var_notes']]

#### Convert new lines to line breaks
#for i in range(len(a['Text_Supplement'])):
#    a['Text_Supplement'][i] = a['Text_Supplement'][i].replace('\n','<br />\n')
#for var_type in v:
#    for var_dict in v[var_type]:
#        for i in range(len(var_dict['Var_Notes'])):
#            var_dict['Var_Notes'][i] = var_dict['Var_Notes'][i].replace('\n','<br />\n')
            
#### If Var_Notes are too long, truncate it so the pdf can be created.
for var_type in v:
    for var in v[var_type]:
        nparag = len(var['var_notes'])
        if nparag > PARAG_LIMIT:
            print 'WARNING: Truncating Var_Notes for %s' % (var['name'])
            var['var_notes'] = var['var_notes'][:PARAG_LIMIT]
            var['var_notes'].append('NOTE: Var_Notes truncated. See NC file for full description.')
            
            
###################################################################################
################################## Write Document #################################
###################################################################################
# Use variables "a" (global attributes), "dims" (dimensions) and "v" (variable 
# attributes, sorted into a dictionary by Var_Type) to populate the pdf

# Create document 
doc = SimpleDocTemplate(fn_out ,pagesize=letter,
                        rightMargin=0.75*inch,leftMargin=0.75*inch,
                        topMargin=72,bottomMargin=18)
Story = []
t = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
styles.add(ParagraphStyle(name='Small', fontSize=8))
styles.add(ParagraphStyle(name='Smallish', fontSize=9))

####################### Introduction ##############################
# Title
text = 'ICON %s' % (a['data_type'].split('>')[-1][1:])
Story.append(Paragraph(text, styles["Title"]))
Story.append(Spacer(1, 12))

# Error notes
if error_notes:
    for text in error_notes:
        Story.append(Paragraph('<font color="red">%s</font>'%text, styles["Normal"]))
        Story.append(Spacer(1, 4))

# Introduction
text = """
This document describes the data product for %s, which is in NetCDF4 format.
""" % (a['description'])
Story.append(Paragraph(text, styles["Justify"]))
Story.append(Spacer(1, 12))

# Text Supplement: one paragraph per string
for text in a['text_supplement']:
    Story.append(Paragraph(text, styles["Justify"]))
    Story.append(Spacer(1, 12))



# Introduction
text = """NetCDF files contain <b>variables</b> and the <b>dimensions</b> over which 
those variables are defined. First, the dimensions are defined, then all variables in the 
file are described."""
Story.append(Paragraph(text, styles["Justify"]))
Story.append(Spacer(1, 12))

####################### Dimensions ##############################
text1 = 'Dimensions'
Story.append(Paragraph(text1, styles["Heading1"]))
Story.append(Spacer(1,6))

text = """
The dimensions used by the variables in this file are given below, along with nominal sizes. Note that the size 
may vary from file to file. For example, the "Epoch" dimension, which describes the number of time samples 
contained in this file, will likely have a varying size.
"""
Story.append(Paragraph(text, styles["Justify"]))
Story.append(Spacer(1, 12))

# Initialize table with row headers
tab_data = []
tab_data.append([Paragraph('<b>Dimension Name</b>',styles['Normal']),
                 Paragraph('<b>Nominal Size</b>',styles['Normal']),
                ])
for dim in dims.keys(): # Construct this row of the table.
    # Build table entries as a Paragraph.
    name_p = Paragraph('<font face="Courier">%s</font>' % (dim), styles["Smallish"])
    dims_p = Paragraph('%s'% (dims[dim]), styles["Smallish"])
    row = [name_p, dims_p]
    tab_data.append(row)           

tab = Table(tab_data, colWidths=[4*inch, 1.1*inch],repeatRows=1)
tab.setStyle(TableStyle([('FONT',(0,1),(0,-1),"Courier"), # filename in Courier
                         ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black), # Grid
                         ('BOX', (0,0), (-1,-1), 0.25, colors.black), # Grid
                         ('VALIGN', (0,0), (-1,-1), "TOP"), # Grid
                        ]))
Story.append(tab)
Story.append(Spacer(1, 12))



####################### Variables ##############################
Story.append(PageBreak()) # new page
text1 = 'Variables'
Story.append(Paragraph(text1, styles["Heading1"]))
Story.append(Spacer(1, 6))

text = """
Variables in this file are listed below. First, "data" variables are described, 
followed by the "support_data" variables, and finally the "metadata" variables. The variables classified as "ignore_data" 
are not shown.
"""
Story.append(Paragraph(text, styles["Justify"]))
Story.append(Spacer(1, 6))

# Make sure variable types are in order: 'data', 'support_data', 'metadata', then any other ones
var_type_ordered = [x for x in ['data','support_data','metadata'] if x in v.keys()]
# Add any other ones (are we expected anything other than those 3?)
for var_type in v.keys():
    if var_type not in var_type_ordered:
        var_type_ordered.append(var_type)
        
for var_type in var_type_ordered:
    text = var_type
    Story.append(Paragraph(text, styles["Heading2"]))
    Story.append(Spacer(1, 6))
    
    # Initialize table with row headers
    tab_data = []
    tab_data.append([Paragraph('<b>Variable Name</b>',styles['Normal']),
                     Paragraph('<b>Description</b>',styles['Normal']),
                     Paragraph('<b>Units</b>',styles['Normal']),
                     Paragraph('<b>Dimensions</b>',styles['Normal']),
                    ])
    for var in v[var_type]: # Construct this row of the table.
        # Build table entries as a Paragraph.
        # First, build the multiple paragraph entry in "Description" cell
        desc_col = [Paragraph(var['catdesc'], styles["Smallish"])] # Start with short description
        ## HACK
        desc_col.append(Spacer(1,6))
        for text in var['var_notes']: # Append paragraphs for each string in Var_Notes
           desc_col.append(Spacer(1,6))
           desc_col.append(Paragraph(text, styles["Small"]))
        # Second, create other cells
        name_p = Paragraph('<font face="Courier">%s</font>' % (var['name']), styles["Smallish"])
        units_p = Paragraph(var['units'], styles["Smallish"])
        dims_s = ', '.join(var['dims'])
        dims_p = Paragraph('<font face="Courier">%s</font>'%dims_s , styles["Smallish"])
        # Put all cells together in a row and add to the table
        row = [name_p, 
               desc_col, 
               units_p, 
               dims_p]
        tab_data.append(row)           

    tab = Table(tab_data, colWidths=[1.8*inch,3.1*inch,0.6*inch,1.3*inch],
                repeatRows=1)
    tab.setStyle(TableStyle([('INNERGRID', (0,0), (-1,-1), 0.25, colors.black), # Grid
                             ('BOX', (0,0), (-1,-1), 0.25, colors.black), # Grid
                             ('VALIGN', (0,0), (-1,-1), "TOP"), # Grid
                            ]))
    Story.append(tab)
    Story.append(Spacer(1, 20))
    

                     
###################### Final Stuff #####################
# Generation notes
text1 = 'This document was automatically generated on <font face="Courier">%s</font> using the file:'%t
text2 = '<font face="Courier">%s</font>' % (fn_in.split('/')[-1])
text3 = 'Software version: <font face="Courier">%s</font>' % (a['software_version'])
Story.append(Paragraph(text1, styles["Small"]))
Story.append(Paragraph(text2, styles["Small"]))
Story.append(Paragraph(text3, styles["Small"]))
Story.append(Spacer(1, 12))

doc.title = a['description']
doc.build(Story)

print 'Created %s' % fn_out