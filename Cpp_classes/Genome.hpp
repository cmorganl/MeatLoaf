#ifndef GENOME_H
#define GENOME_H

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

class Genome {
    protected:
        int longest_contig;
        Genome();
//        ~Genome();
        virtual int record_header( string, int , int* ) =0;
        virtual int record_sequence( string, int ) =0; 
    public:
        int find_longest_contig();
        char **header;
        char **sequence;
        int N_contigs;
        int genome_length;
        char ** retrieve_headers() { return header; };
        char ** retrieve_sequences() { return sequence; };
};

#endif
