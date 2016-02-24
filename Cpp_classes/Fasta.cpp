/*
 * Fasta.cpp (c) 2015 Connor Morgan-Lang <c.morganlang@gmail.com>
 * Hallam Lab, Department of Microbiology and Immunulogy, UBC
 * ------------------------------------------
 * Last modified: 23 February 2016 (CML)
 * ------------------------------------------
 */

#include "Fasta.hpp"

Fasta::Fasta(std::string input) {
    sequences = (char **) malloc (sizeof(char *));
    longest_contig = 0;
    genome_length = 0;
    N_contigs = 0;
    header_base.seq_length.push_back(0);

    fasta_file = new ifstream(input.c_str());

    if ( !fasta_file->is_open() ) {
        cerr << "Unable to open '" << input << "' for reading. Exiting now!" << endl;
        exit(0);
    }
}

Fasta::Fasta(void) { clear(); }

void Fasta::clear(void) {
    int k;
    for (k = 0; k < N_contigs; k++) {
        header_base.name.pop_back();
        header_base.seq_length[k] = 0;
        sequences[k] = NULL;
    }
    longest_contig = 0;
    genome_length = 0;
    N_contigs = 0;
}

Fasta::~Fasta() {
    fasta_file->close();
    int k;
    for (k = 0; k < N_contigs; k++) {
        free(sequences[k]);
    }
    free(this->sequences);

}

// copy constructor:
Fasta::Fasta( const Fasta& other) : header_base(other.header_base), sequences(other.sequences) {
    std::cout << "Copy constructor activated" << std::endl;
}

int Fasta::find_longest_contig() {
    int x;
    int size = 0;
    cout << "Number of contigs = " << N_contigs << endl;
    std::cout << "Finding the longest contig... ";
    for (x = 0; x < N_contigs; x++) {
        if ( sequences[x] == NULL) {
            fprintf(stderr, "ERROR: A sequence was improperly loaded into FASTA!\n");
            cout << "Index " << x << " (" << header_base.name[x] << ") is NULL." << endl;
            exit(-1);
        }
        size = strlen(sequences[x]);
        if (size != header_base.seq_length[x] )
            fprintf(stderr, "ERROR: Sequences and header_base are not matching!\n");

        if (size > longest_contig)
            longest_contig = size;
    }
    std::cout << "done." << endl;
    return longest_contig;
}

int Fasta::find_sequence_length(int number) {
    if ( number <= N_contigs )
        return this->header_base.seq_length[number];
    else {
        fprintf(stderr, "ERROR: header_base.seq_length index out of range!\n");
        exit(-1);
    }
}

int Fasta::record_header( string line ) {
    //TODO: Check to make sure there are no duplicate headers
    header_base.name.push_back(line);
    header_base.seq_length.push_back(0);
    return 0;
}

int Fasta::record_sequence( string line, int line_len) {
    genome_length += line_len - 2;
    header_base.seq_length[N_contigs] += line_len - 2;
    sequences = (char **) realloc (sequences, (genome_length + N_contigs) * sizeof(char *));
    if ( sequences[N_contigs] == NULL) {
        sequences[N_contigs] = (char *) malloc (line_len * sizeof(char));
        strcpy(sequences[N_contigs], line.c_str());
    }
    else {
        sequences[N_contigs] = (char *) realloc (sequences[N_contigs], (strlen(sequences[N_contigs]) + line_len) * sizeof(char));
        strcat(sequences[N_contigs], line.c_str());
    }
    return 0;
}

int Fasta::parse_fasta(int min_length) {
    string line;
    int line_len;
    std::string header;
    sequences[N_contigs] = NULL;

    if (fasta_file->good()) {
        getline( *fasta_file, line);
    }
    else {
        fprintf(stderr, "The fasta file cannot be parsed!\n");
        exit(0);
    }
    while ( fasta_file->good() ) {
        if (line.length() == 0)
            ;
        else {
            line_len = line.length() + 2;
            if (line.at(0) == '>') {
                if ( find_sequence_length(N_contigs) >= min_length ) {
                    record_header(header);
                    N_contigs++;
                }
                else {
                    header_base.seq_length[N_contigs] = 0;
                    sequences[N_contigs] = NULL;
                }
                header = line;
            }
            else
                record_sequence( line, line_len );
        }
        getline( *fasta_file, line );
    }
    if ( find_sequence_length(N_contigs) >= min_length ) {
        record_header(header);
        N_contigs++;
    }
    else
        sequences[N_contigs] = NULL;

    if ( (signed)header_base.name.size() != N_contigs )
        return 1;
    else
        return 0;
}
