#include "Genome.hpp"

Genome::Genome() {
    header = (char **) malloc (sizeof(char *));
    sequence = (char **) malloc (sizeof(char *));
    longest_contig = 0;
    genome_length = 0;
    N_contigs = 0;
}
/*
Genome::~Genome () {
    free(header);
    free(sequence);
}
*/
int Genome::find_longest_contig() {
    int x;
    int size = 0;
    for (x = 0; x < N_contigs; x++) {
        size = strlen(sequence[x]);
        if (size > longest_contig)
            longest_contig = size;
    }
    return longest_contig;
}
