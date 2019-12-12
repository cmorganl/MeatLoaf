/*
 * Fasta.cpp Connor Morgan-Lang <c.morganlang@gmail.com>
 * Hallam Lab, Department of Microbiology and Immunulogy, UBC
 * ------------------------------------------
 * Last modified: April 11 2016 (CML)
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

int Fasta::writeNx(string output, bool verbose) {
    /*
     * A function for writing a csv for generating Nx curves
     * Columns are Genome_proportion, Contig_size
     */
    if (verbose)
        cerr << "Sorting sequence lengths... ";
    std::sort (header_base.seq_length.begin(), header_base.seq_length.end());
    if (verbose)
        cerr << "done." << endl;

    float x = 0.0;
    signed long prop_counter = N_contigs;
    signed long seq_mass = header_base.seq_length[prop_counter];
    if (verbose)
        cerr << "Sequence mass initialized to longest sequence: " << seq_mass << "bp" << endl;

    std::ofstream Nx_out;
    if (output.size() > 1) {
        if (verbose)
            cerr << "Opening output \"" << output << "\"" << endl;
        Nx_out.open(output.c_str());
        Nx_out << "Genome_proportion" << "," << "Contig_size" << endl;
    }
    else {
        cerr << "No output provided - printing to stdout" << endl;
        cerr << endl;
        cout << "Genome_proportion" << "," << "Contig_size" << endl;
    }

    while ( x <= 100 ) {
        while (seq_mass < ((x / 100) * genome_length) && prop_counter > 0) {
            prop_counter--;
            if (header_base.seq_length[prop_counter] == 0) {
                cerr << "Sequences with 0 length were recorded at index " << prop_counter << endl;
                exit(3);
            }
            seq_mass += header_base.seq_length[prop_counter];

        }
        if (output.size() > 1) {
            Nx_out << x / 100 << "," << header_base.seq_length[prop_counter] << endl;
        }
        else
            cout << x / 100 << "," << header_base.seq_length[prop_counter] << endl;
        x++;
    }
    if (output.size() > 1) {
        Nx_out.close();
    }
    return 1;
}

long Fasta::find_longest_contig() {
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

long Fasta::find_sequence_length(long number) {
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

int Fasta::record_sequence( string line ) {
    int line_len = line.length();
    header_base.seq_length[N_contigs] += line_len;
    sequences = (char **) realloc (sequences, (genome_length + line_len) * sizeof(char *));
    if ( sequences[N_contigs] == NULL ) {
        sequences[N_contigs] = (char *) malloc ( (line_len+2) * sizeof(char));
        strcpy(sequences[N_contigs], line.c_str());
        sequences[N_contigs + 1] = NULL;
    }
    else {
        sequences[N_contigs] = (char *) realloc (sequences[N_contigs],
                                                 (strlen(sequences[N_contigs]) + (line_len+2)) * sizeof(char));
        strcat(sequences[N_contigs], line.c_str());
    }
    return 0;
}

int Fasta::parse_fasta(int min_length) {
    string line;
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
            if (line.at(0) == '>') {
                if ( find_sequence_length(N_contigs) >= min_length ) {
                    record_header(header);
                    genome_length += header_base.seq_length[N_contigs];
                    N_contigs++;
                }
                else {
                    header_base.seq_length[N_contigs] = 0;
                    sequences[N_contigs] = NULL;
                }
                header = line;
            }
            else
                record_sequence(line);
        }
        getline( *fasta_file, line );
    }
    genome_length += header_base.seq_length[N_contigs];
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
