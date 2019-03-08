# Motif-mark
Python script to visualize motifs on sequences. One image will be created for all reads.
__Default output file is an SVG image.__ Motif colors are chosen randomly, if you're 
unhappy with colors or they are too similar, simply run the program until you obtain 
desirable colors.

__Usage__ 

python motifMark.py [Options]

__Input flags__

_Required_

-f [Fasta file] - supports multiple fasta sequences

-m [File containing interested motifs] - supports multiple motifs

_Optional_

-p [pdf] - Output drawn image as a pdf file

-n [png] - Output drawn image as a png file

__Example__

python motifMark.py -f Figure_1.fasta -m Fig_1_motifs.txt

__Coming soon__

-c [Color palette] - file containing RGB or hexadecimal colors for motifs


