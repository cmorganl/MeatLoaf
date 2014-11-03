#ifndef FASTA_H
#define FASTA_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include "Genome.hpp"

using namespace std;

class FastaParser : public Genome {
    public:
        ifstream *fasta_file;
        int header_Nchar;
        FastaParser( char * );
        int parse_fasta();
        int record_header( string, int, int * );
        int record_sequence( string, int );
};

#endif
