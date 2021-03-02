#!/usr/bin/env python

#####################################
##########     SET-UPS     ##########
#####################################

import argparse
import itertools
import re
import cairo
import seaborn as sns 


def get_args():
    """
    Will set a command line option to run arguments.
    """
    # Parser main
    parser = argparse.ArgumentParser(description = "Creates a gene map image of motif binding sites for \
        a given sequence or set of sequences. One SVG image is created for each input <file.fasta>. \
        The <file.fasta> name will be used to name the output <image.svg>. Image will be saved \
        in the directory given.")
    
    # Parser arguments
    parser.add_argument('-f','--FASTAfile', help = "Absolute/path/to/<file.fasta>.", required = True)
    parser.add_argument('-m','--MOTIFfile', help = "Absolute/path/to/<motif.txt>.", required = True)

    return parser.parse_args()
args = get_args()

# Get filenames and directory
FASTA_file = args.FASTAfile
MOTIF_file = args.MOTIFfile

# IUPAC dictionary of ambiguous base changes
IUPAC_dict = {
    # Use lowercase letters
    "a":["a"            ],
    "c":[    "c"        ],
    "g":[        "g"    ],
    "t":[           "u","t"],
    "u":[           "u","t"],
    "w":["a",       "t","u"],
    "s":[    "c","g"    ],
    "m":["a","c"        ],
    "k":[       "g","t","u"],
    "r":["a",    "g",   ],
    "y":[    "c",   "t","u"],
    "b":[   "c","g","t","u"],
    "d":["a",   "g","t","u"],
    "h":["a","c",   "t","u"],
    "v":["a","c","g",   ],
    "n":["a","c","g","t","u"],
    "z":[               ],
}

### Image presets
GENE_GROUP_HEIGHT = 100 # The space for which a gene group fits into
X0,Y0 = 10,10 # Coordinate offsets


#####################################
##########     CLASSES     ##########
#####################################


class GeneGroup:
    # Data
    def __init__(self, gene_num, header, seq, seq_length):
        self.gene_num = gene_num
        self.header = header
        self.seq = seq
        self.seq_length = seq_length


class FastaHeader:
    # Data
    def __init__(self, gene_num, header):
        self.gene_num = gene_num
        self.header = header
    
    # Behavior
    def draw(self, context, header, gene_num):
        """
        Draws fasta header on image.
        """
        y = Y0 + 100 * gene_num
        context.set_source_rgb(0,0,0) # black text
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(12)
        context.move_to(X0,y)
        context.show_text(header)

class Gene:
    # Data
    def __init__(self, gene_num, seq_length):
        self.gene_num = gene_num
        self.seq_length = seq_length
    
    # Behavior
    def draw(self, context, gene_num, seq_length):
        """
        Draws entire gene (i.e., intronic regions)
        """
        height = 5
        #x = X0 + seq_length
        y = Y0 + float(GENE_GROUP_HEIGHT / 2.5) + (GENE_GROUP_HEIGHT * gene_num)
        context.set_source_rgb(0,0,0) # black color
        cr.rectangle(X0, y, seq_length, height)
        cr.fill()
        

class Exon:
    # Data
    def __init__(self, gene_num, exon_start, exon_length):
        self.gene_num = gene_num
        self.exon_start = exon_start
        self.exon_length = exon_length

    # Behavior
    def draw(self, context, gene_num, exon_start, exon_length):
        """
        Draws exons on gene.
        """
        height = 35
        x1 = X0 + exon_start
        y = (GENE_GROUP_HEIGHT // 2.5) + (GENE_GROUP_HEIGHT * gene_num) - 5
        context.set_source_rgb(0,0,0) # black color
        cr.rectangle(x1, y, exon_length, height)
        cr.fill()
        

class Motif:
    # Data
    def __init__(self, name, length, color, combos, motif_coords):
        self.name = name
        self.length = length
        self.color = color
        self.combos = combos
        self.motif_coords = motif_coords
        
    # Behavior
    def draw(self, context, gene_num, motif_pos, length, color):
        """
        Draws motifs (by color) according to its sequence length.
        """
        height = 35
        x1 = X0 + motif_pos
        y2 = (GENE_GROUP_HEIGHT // 2.5) + (GENE_GROUP_HEIGHT * gene_num) - 5 - GENE_GROUP_HEIGHT
        context.set_source_rgba(color[0], color[1], color[2], 0.8)
        context.rectangle(x1, y2, length, height)
        context.fill()
        

class Legend:
    # Data
    def __init__(self, motif_num, name, length, color):
        self.motif_num = motif_num
        self.name = name
        self.length = length
        self.color = color
    
    # Behavior
    def label_legend(self, context, name, motif_num):
        x = (X0 + 100) * motif_num
        y = IMG_HEIGHT - 45//2
        context.set_source_rgb(0,0,0) # black text
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(12)
        context.move_to(x,y)
        context.show_text(name)

    def color_legend(self, context, motif_num, length, color):
        """
        Draws legend at the bottom of the image.
        """
        height = 35
        x1 = ( (X0 + 100) * motif_num ) - 15
        y = IMG_HEIGHT - 45
        context.set_source_rgba(color[0], color[1], color[2], 0.8)
        context.rectangle(x1, y, length, height)
        context.fill()


#####################################
##########   FASTA STUFF   ##########
#####################################


### Read in file_fasta to save objects as GeneGroup class
GENEgroup_obj_list = []
with open(FASTA_file, 'r') as fasta_fh:
    gene_num = 0
    LN = 0
    for line in fasta_fh:
        LN += 1

        # Make headers GeneGroup object
        if line.startswith(">") == True:
            seq_list = []
            gene_num += 1
            header = line.strip()
            obj = GeneGroup(gene_num, header, "", "")
            
            # Save object to list
            GENEgroup_obj_list.append(obj)
            continue
            
        # Concatenate sequence lines
        seq = line.strip()
        seq_list.append(seq)
        long_seq = "".join(seq_list)

        # Get length of sequence lin
        seq_length = len(long_seq)

        # Save to GeneGroup object
        obj.seq = long_seq
        obj.seq_length = seq_length

### Get new objects from GeneGroup objects
FASTAheader_obj_list = []
GENE_obj_list = []
EXON_obj_list = []
for count, obj in enumerate(GENEgroup_obj_list):
    
    # Add FastaHeader objects to list
    obj1 = FastaHeader(count, obj.header)
    FASTAheader_obj_list.append(obj1)

    # Add Gene objects to list
    obj2 = Gene(count, obj.seq_length)
    GENE_obj_list.append(obj2)

    # Use regex to find exons (all uppercase letters)
    exon = re.finditer(r"[A-Z]", obj.seq)

    # Get exon start position and length
    exon_pos = []
    for region in exon:
        pos = region.start() + 1 # adjust to make 1-based
        exon_pos.append(pos)
    
    # Add Exon objects to list
    obj3 = Exon(count, exon_pos[0], len(exon_pos))
    EXON_obj_list.append(obj3)


#####################################
##########   MOTIF STUFF   ##########
#####################################


### Read in motif.txt and add to motif object list
MOTIF_obj_list = []
with open(MOTIF_file, 'r') as motif_fh:
    LN = 0

    # Loop through file to get motif sequences
    for line in motif_fh:
        LN += 1
        motif_name = line.strip().lower()

        # Save as motif object
        motif_obj = Motif(motif_name, len(motif_name), "", "", "")

        # Add motif object to list
        MOTIF_obj_list.append(motif_obj)

### Generate color palette
palette = sns.color_palette("Set2", len(MOTIF_obj_list))

### Get new info for motif objects
for count, motif_obj in enumerate(MOTIF_obj_list):
    # Assign colors for each motif_obj
    motif_obj.color = palette[count]

    # Create lambda function called groups
    # Makes an intertools object for each motif by substituting characters 
    groups = itertools.groupby(motif_obj.name, lambda char:char not in IUPAC_dict)
    splits = []

    # Loop through each itertools group object
    for b, group in groups:

        # Motif contains ambiguous character
        if b:
            # Add translated base to end of sequence
            splits.extend([[g] for g in group])
        
        # No ambiguous character in motif
        else:

            # Add nucleotide base to end of sequence
            for nuc in group:
                splits.append(IUPAC_dict[nuc])

    # Create dictionary and add translated sequences as values
    combos = ["".join(p) for p in itertools.product(*splits)]

    # Add to combos to motif object
    motif_obj.combos = combos

    ### Get all motif coordinates from each gene and save to motif_obj
    # Loop through each gene object to save motif coordinates to motif_obj_dict
    motif_pos_dict = {}
    for count2, gene_obj in enumerate(GENEgroup_obj_list):
        sequence = gene_obj.seq

        # Loop through each motif combination to search in sequence
        motif_pos = []
        for motif_combo in motif_obj.combos:
            pattern = str("(?:" + motif_combo + ")")

            # Use regex to find motif coordinate in current sequence
            coordinate = re.finditer(pattern, sequence)

            # Get starting positions of motifs
            for pos in coordinate:
                pos = pos.start() + 1 # adjust to make 1-based
                motif_pos.append(pos)
        
        # Add coordinate start positions to dictionary
        # key = gene_num; value = list of motif start positions
        motif_pos_dict[count2 + 1] = motif_pos

    # Save motif coordinates to motif_obj
    motif_obj.motif_coords = motif_pos_dict


#####################################
##########     LEGEND      ##########
#####################################


### Loop through motif objects to get Legend objects
LEGEND_obj_list = []
for count, motif_obj in enumerate(MOTIF_obj_list):
    obj = Legend(count + 1, motif_obj.name, motif_obj.length, motif_obj.color)
    LEGEND_obj_list.append(obj)


#####################################
##########    IMG STUFF    ##########
#####################################


def name_image(file_name):
    """
    Names image.svg from path given to args. 
    Example:

    args.FASTAfile = "full/path/to/the/file_name.fasta"

    name_image(args.FASTAfile) = "file_name.svg"
    """
    # Make image name the same as filename
    image_name = "%s" % FASTA_file
    image_name = image_name.split("/")
    image_name = image_name[-1]
    image_name = image_name.split(".")
    image_name = str(image_name[0]) + ".svg"
    return image_name
IMG_NAME = name_image(FASTA_file)

def get_dims(GENEgroup_obj_list):
    """
    Sets image dimensions.
    """
    lengths = []
    for genegroup_obj in GENEgroup_obj_list:
        lengths.append(genegroup_obj.seq_length)
        
    WIDTH = 30 + max(lengths) * 1.25
    HEIGHT = len(lengths) * 100 + 50 # Add extra space for legend
    return WIDTH, HEIGHT
IMG_WIDTH, IMG_HEIGHT = get_dims(GENEgroup_obj_list)

### Begin drawing image
with cairo.SVGSurface(IMG_NAME, IMG_WIDTH, IMG_HEIGHT) as surface:
    # Start object for SVG surface
    cr = cairo.Context(surface)

    # Draw FastaHeaders
    for fasta_obj in FASTAheader_obj_list:
        fasta_obj.draw(cr, fasta_obj.header, fasta_obj.gene_num)
    
    # Draw Genes
    for gene_obj in GENE_obj_list:
        gene_obj.draw(cr, gene_obj.gene_num, gene_obj.seq_length)
    
    # Draw Exons
    for exon_obj in EXON_obj_list:
        exon_obj.draw(cr, exon_obj.gene_num, exon_obj.exon_start, exon_obj.exon_length)
    
    # Draw Motifs
    for motif_obj in MOTIF_obj_list:
        color = motif_obj.color
        motif_coord_dict = motif_obj.motif_coords

        # Loop through motif coordinates dictionary
        for item in motif_coord_dict:
            gene_num = item
            coordinates = motif_coord_dict[item]

            # Loop through each motif coordinate to draw on image
            for coordinate in coordinates:
                motif_obj.draw(cr, gene_num, coordinate, motif_obj.length, color)
        
    # Draw Legend
    for legend_obj in LEGEND_obj_list:
        legend_obj.label_legend(cr, legend_obj.name, legend_obj.motif_num)
        legend_obj.color_legend(cr, legend_obj.motif_num, legend_obj.length, legend_obj.color)
        