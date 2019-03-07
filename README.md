# Motif-mark
Python script to visualize motifs on sequences. An image will be created for each read.
Default output file is an SVG image. Motif colors are chosen randomly, if you're 
unhappy with colors or they are too similar, simply run the program until you obtain 
desirable colors.

__Usage__ 
python motifMark.py [options]

__Input flags__

_Required_

-f [Fasta file] - supports multiple fasta sequences

-m [file containing interested motifs] - supports multiple motifs

__Coming Soon__
-c [color palette] - file containing RGB or hexadecimal colors for motifs

-p [Output as png] 

-j [Output as jpeg]

