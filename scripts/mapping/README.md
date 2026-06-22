Read mapping and mapping-statistics generation using the vg toolkit

For each reference substrate, reads were mapped using the vg toolkit. Linear references were first indexed using vg autoindex, with separate workflows for short-read and long-read mapping.

For short-read mapping, indexes were generated as follows:

vg autoindex \
    --workflow giraffe \
    -r FASTA.fasta \
    -p PREFIX

For long-read mapping, indexes were generated using the lr-giraffe workflow:

vg autoindex \
    --workflow lr-giraffe \
    -r FASTA.fasta \
    -p PREFIX

Short-read paired-end data were mapped using vg giraffe:

vg giraffe \
    --hard-hit-cap 9000 \
    -p \
    -t 32 \
    -Z PREFIX.gbz \
    -d PREFIX.dist \
    -m PREFIX.shortread.withzip.min \
    -f R1.fastq.gz \
    -f R2.fastq.gz \
    > Alignments.gam

PacBio HiFi long reads were mapped using vg giraffe with the HiFi preset:

vg giraffe \
    -b hifi \
    -p \
    -t 32 \
    -Z PREFIX.gbz \
    -d PREFIX.dist \
    -m PREFIX.longread.withzip.min \
    -z PREFIX.longread.zipcodes \
    -f READS.fastq.gz \
    > Alignments.gam

Mapping statistics were generated from the resulting GAM files using vg stats, which is part of the vg toolkit:

vg stats -a Alignments.gam

The resulting summaries were used to calculate mapping rates and alignment-quality statistics across the different reference substrates.For mapping to a pangenome graph, the same indexing and mapping workflow was used, except that the linear FASTA input was replaced with a graph input during indexing. Specifically, -r FASTA.fasta was replaced with either -g GFA_file or -G GBZ_file, depending on the graph format used.
