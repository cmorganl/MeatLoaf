#ifndef FASTA_H
#define FASTA_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include <string>

using namespace std;

class Fasta {
protected:
    int record_header( std::string, int, int * );
    int record_sequence( std::string, int );
    //vector Fastq2Fasta( char * ); TODO:CONNOR needs to be implemented
public:
    Fasta( std::string );
    Fasta( void );
    ~Fasta();
    Fasta( const Fasta& other);
    void clear( void );

    char **header;
    char **sequence;
    int N_contigs;
    int genome_length;
    int longest_contig;
    ifstream *fasta_file;
    int header_Nchar;

    char ** retrieve_headers() { return header; };
    char ** retrieve_sequences() { return sequence; };
    int parse_fasta();
    int find_longest_contig();
};

#endif