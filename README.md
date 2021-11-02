# airen_r
airen: A R package to optimize the performance of Restriction digestion to improve the resolution of WGS


# Requirements:

- R environment (>4.0).
- Reference genome in FASTA file (.fa)
- Genomic aligner (like HISAT2).

# Functions and options

In general, the main package is generally adressed to define the best combination of restriction enzymes over a genome in order to increase the resolution of WGS sequencing. The functions of the package are defined by:

- restrictionSimulation: Performance of in-silico restriction enzymes using SimRAD package using our function restrictionSimulation (you have to stablish the range of desired size of fragments, and multiple options are defined to choose your action: 'finding', 'combination' and 'replicate'.

'finding': Use a database of 400 enzyme to predict the best enzyme to obtain the most quantity of fragments in the desired size range from the reference genome.
'combination': Use two chosen enzymes (presents in the airen enzyme database) to perform the in-silico digestion.
'replicate': Use the chosen enzyme(s) to repeat the same in-silico digest in order to predict accurately the number of expected fragments.

- randomGenomeFragmentation: Performance of random genome fragmentation using the reference genome and stablishing the number of repetitions to have a good representation of control data.
- calculatePositionFragment: Function that allows to locate the aligned fragments in the reference genome, in order to identify the cut sites.

# Workflow

The main workflow is defined by that steps:

restrictionSimulation function

- Defining the best enzyme to increase the number of expected framgents from a desired range using restrictionSimulation and 'finding' condition.
- Performing the combination of the two best enzymes with the most quantity of desired fragments to show the potential of the combination using 'combination' condition.
- Repeat n times the restriction to increase the accuracy of the prediction using 'replicate condition'.

randomGenomeFragmentation function

- Generate n random genome fragmentation in order to generate a control to compare the efficiency of the restriction enzyme.

Alignment with some aligner like HISAT2 and generate a SAM file.

calculatePositionFragment function

- Identify and compare the control random genome fragmentation to the resulted fragments from the in-silico digestion using the SAM file.
