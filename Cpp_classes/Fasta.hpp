#ifndef FASTA_H
#define FASTA_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include <string>
#include <vector>
#include <algorithm>

using namespace std;

struct header {
    std::vector<std::string> name;
    std::vector<int> seq_length;
};

class Fasta {
protected:
    int record_header( std::string );
    int record_sequence( std::string );

public:
    // Initialization functions
    Fasta( std::string );
    Fasta( void );
    ~Fasta();
    Fasta( const Fasta& other);
    void clear( void );

    // Class objects
    header header_base;
    char **sequences;
    long int N_contigs;
    long int genome_length;
    int longest_contig;
    ifstream *fasta_file;

    // Class functions
    std::vector<std::string> retrieve_headers() { return header_base.name; };
    char ** retrieve_sequences() { return sequences; };
    int parse_fasta(int min_length);
    int find_longest_contig();
    int find_sequence_length(int number);
    int writeNx(std::string output, bool verbose);
};

#endif