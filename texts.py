GRAPH_DESCRIPTION = """
The heatmap below  displays  the comparison of the alignments results.
Every selected aligner is displayed on the left in the heatmap and all
proteins from the input network are listed as columns. The agreement
level between the aligned proteins is represented through a color
palette that ranges from light green to dark green:
*  Dark green indicates a total agreement, meaning that all aligners have matched the protein to
the same protein in the output network.
* Light green color means that all aligners matched the corresponding protein to a different one in the
output network. Thus, a vertical column in dark green shows a consensus
among all aligners.
* The color corresponding to unaligned proteins is
gray.

In addition, when pointing the cursor in the columns on the
heatmap a window with the following information is shown: 
* *Origin*: the protein name in the input network.
* *Target*: the protein name in the output network to whom the origin protein has been assigned by the aligner.
* *Aligner*: the corresponding aligner.
"""

TABLE_DESCRIPTION = """
In the page size below the number of proteins (columns) to visualize in
the heatmap can be selected.  Also, the table displays the alignment
results. For every row, the first column is the protein name in the
input network. Its corresponding matching for every aligner are listed
in the other columns. Last column provides a consensus score which is
the mean of agreement among all aligners.
"""