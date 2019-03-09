#!/usr/bin/env python3
import cairo
import math
import argparse
import re
import random

def get_arguments():
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Python script to visualize motifs on sequences. One image will be created for all reads. Default output file is an SVG image. Motif colors are chosen randomly, if you\'re unhappy with colors or they are too similar, simply run the program until you obtain desirable colors.')
        parser.add_argument("-f", "--file", help="Input fasta file", type = str, required = True)
        parser.add_argument("-m", "--motif", help="Input motif file", type = str, required = True)
        parser.add_argument("-p", "--pdf", help="Output image as a PDF file", required = False, action="store_true", default = False)
        parser.add_argument("-n", "--png", help="Output image as a png file", required = False, action="store_true", default = False)
        return parser.parse_args()

def parseMotif(motifFile):
    """ This function takes a file of motifs and returns all of the motifs as a list"""
    with open(motifFile) as motifFile:
        motifs = []
        for line in motifFile:
            line = line.strip().split()
            for motif in line:
                motif = motif.lower()
                motifs.append(motif)
        return(motifs)
    
def parseFasta(fastaFile):
    """ This function takes a fasta file as input and returns the header of each sequence as a key and the fasta sequence
    as the value"""
    with open(fastaFile) as fastFile:
        data = ''
        #Holds header as key and fasta sequence as value
        Seqs = {}
        for line in fastFile:
            line = line.strip()
            if line.startswith(">"):
                if data:
                    Seqs[inKey] = data
                    data = ''
                inKey = line
            elif not line.startswith(">"):
                data += line
            Seqs[inKey] = data 
    return(Seqs)
    
def parseSeqDict(fastaSeq, MotifList):
    """ This function takes a fasta sequence and List of motifs and returns each motif's relative position along with 
    the position of all Exons in the fasta sequence"""
    #Holds relative positions of Exons
    exonPos = []
    #Holds start and end position on line of each motif match as value and original motif 
    motifPos = {}
    exons = re.compile("[ATGCYRSWKMBDHVN]+")
    for motif in MotifList:
        oldMotif = motif
        #IUPAC nucleotide codes being replaced in motifs
        if 'y' in motif:
            motif = motif.replace('y', '[ct]')
        if 'r' in motif:
            motif = motif.replace('r', '[ag]')
        if 's' in motif:
            motif = motif.replace('s', '[gc]')
        if 'w' in motif:
            motif = motif.replace('w', '[at]')
        if 'k' in motif:
            motif = motif.replace('k', '[gt]')
        if 'm' in motif:
            motif = motif.replace('m', '[ac]')
        if 'b' in motif:
            motif = motif.replace('b', '[cgt]')
        if 'd' in motif:
            motif = motif.replace('d', '[agt]')
        if 'h' in motif:
            motif = motif.replace('h', '[act]')
        if 'v' in motif:
            motif = motif.replace('v', '[acg]')
        if 'n' in motif:
            motif = motif.replace('n', '[acgt]')
        if 'u' in motif:
            motif = motif.replace('u', '[t]')
        motifMatch = re.compile(motif, re.IGNORECASE)
        motifPos[oldMotif] = []
        #Find motif sequence iteratively in fasta sequence and append the span as tuple to motif dictionary
        for m in re.finditer(motifMatch, fastaSeq):
            motifPos[oldMotif].append(m.span())
        #Find Exons iteratively in fasta sequence and append the span as tuple to exon list
        for e in re.finditer(exons, fastaSeq):
            exonPos.append(e.span())
    return(exonPos, motifPos)
    
def floatRgb(mag, cmin, cmax):
    """ Return a tuple of floats between 0 and 1 for R, G, and B. """
    # Normalize to 0-1
    try: x = float(mag-cmin)/(cmax-cmin)
    except ZeroDivisionError: x = 0.5 # cmax == cmin
    blue  = min((max((4*(0.75-x), 0.)), 1.))
    red   = min((max((4*(x-0.25), 0.)), 1.))
    green = min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
    return(red, green, blue)
    
args = get_arguments()
file = args.file
file2 = args.motif
# file = "Figure_1.fasta"
# file2 = "Fig_1_motifs.txt"



motifList = parseMotif(file2)
Sequences = parseFasta(file)
width, height = 800, 400 * len(Sequences)
#create the coordinates to display your graphic, desginate output
if args.pdf == True:
    surface = cairo.PDFSurface("output2.pdf",width, height)
elif args.png:
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
else:
    surface = cairo.SVGSurface("output2.svg",width, height)
#create the coordinates you will be drawing on (like a transparency) - you can create a transformation matrix
context = cairo.Context(surface)
#context.scale(width,height) #will set your drawing surface to a 0.0-1.0 scale

previousRGB = []
for motif in motifList:
    #Pick random integers betweeen 0-1 for RGB values
    mag, mini, maxi = random.randint(1,49), random.randint(0,50), random.randint(0,50)
    red, green, blue = floatRgb(mag, mini, maxi)
    #Check if color has been used for a previous motif
    while (red, green, blue) in previousRGB:
        mag, mini, maxi = random.randint(1,49), random.randint(0,50), random.randint(0,50)
        red, green, blue = floatRgb(mag, mini, maxi)
        # print("Picking new RGB values")
    previousRGB.append((red, green ,blue))
    
readNum = 0
for key in Sequences:
    read = Sequences[key]
    Exons, Motifs = parseSeqDict(read, motifList)
    startPoint = 50
    endPoint = 750
    lengthRNA = endPoint - startPoint
    nucLength = lengthRNA / len(read)
    strandYPos = 75 + readNum

    #Draw entire strand as introns
    context.set_line_width(5)
    context.move_to(startPoint,strandYPos + readNum)        #(x,y)
    context.line_to(endPoint,strandYPos + readNum)
    context.stroke()

    #For loop draws all of the exons in Exon list of tuples
    for exon in Exons:
        context.set_line_width(30)
        context.move_to(startPoint + (exon[0] * nucLength), strandYPos + readNum)        #(x,y)
        context.line_to(startPoint + (exon[1] * nucLength), strandYPos + readNum)
        context.stroke()
    
    #List hold previous rgb values for unique motifs.
    #For loop draws all of the motifs in Motif list of tuples.
    motifColor = 0
    for motif in Motifs:
        value = Motifs[motif]
        context.set_source_rgb(previousRGB[motifColor][0], previousRGB[motifColor][1], previousRGB[motifColor][2])
        motifColor += 1
        for span in value:
            context.set_line_width(len(motif) / 2)
            context.move_to(startPoint + (span[0] * nucLength), strandYPos - 15 + readNum)
            context.line_to(startPoint + (span[0] * nucLength), strandYPos + 15 + readNum)
            context.stroke()

    #Drawing legend
    context.set_source_rgb(0, 0, 0)
    context.set_line_width(1)
    #Top legend line
    context.move_to(endPoint - 200, strandYPos + 40 + readNum)
    context.line_to(endPoint, strandYPos + 40 + readNum)
    context.stroke()
    #Bottom legend line
    context.move_to(endPoint - 200, strandYPos + 65 + (len(Motifs) * 20) + readNum)
    context.line_to(endPoint, strandYPos + 65 + (len(Motifs) * 20) + readNum)
    context.stroke()
    #Right legend line
    context.move_to(endPoint, strandYPos + 40 + readNum)
    context.line_to(endPoint, strandYPos + 65 + (len(Motifs) * 20) + readNum)
    context.stroke()
    #Left legend line
    context.move_to(endPoint - 200, strandYPos + 40 + readNum)
    context.line_to(endPoint - 200, strandYPos + 65 + (len(Motifs) * 20) + readNum)
    context.stroke()
    #Legend Title
    context.set_source_rgb(0, 0, 0)
    context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(14)
    context.move_to(endPoint - 120, strandYPos + 55 + readNum)
    context.show_text("Motifs")

    logo = 0
    #For loop is used to draw logos in legend
    for color in previousRGB:
        context.set_source_rgb(color[0], color[1], color[2])
        context.rectangle(endPoint - 180, (strandYPos + 60) + (logo) + readNum, 13, 13)
        context.fill()
        logo += 20

    motifText = 0
    #For loop is used to write motif text in legend
    for motif in Motifs:
        context.set_source_rgb(0, 0, 0)
        context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(13)
        context.move_to(endPoint - 160, (strandYPos + 68.5) + (motifText) + readNum)
        context.show_text("- " + motif)
        motifText += 20

    context.set_source_rgb(0, 0, 0)
    context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(10)
    context.move_to(startPoint, strandYPos - 25 + readNum)
    context.show_text(key)
    context.move_to(startPoint, strandYPos + 25 + readNum)
    context.show_text("0")
    context.move_to(endPoint - 10, strandYPos + 25 + readNum)
    context.show_text(str(len(read)))
    
    readNum += 150

if args.png:
    surface.write_to_png("output2.png")
else:
    surface.finish()
print("Image created")