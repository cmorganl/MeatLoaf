#include "FastaParser.hpp"

FastaParser::FastaParser( char * input ) {
    header_Nchar = 0;
    fasta_file = new ifstream(input);
    if ( !fasta_file->is_open() ) {
        cerr << "Unable to open " << input << "for reading. Exiting now!" << endl;
        exit(0);
    }
}

int FastaParser::record_header( string line, int line_len, int * header_Nchar ) {
    * header_Nchar += line_len;
    header = (char **) realloc (header, * header_Nchar * sizeof(char));
    header[N_contigs] = (char *) malloc (line_len * sizeof(char));
    strcpy( header[N_contigs], line.c_str());
    sequence[N_contigs] = NULL;
    return 0;
}

int FastaParser::record_sequence( string line, int line_len) {
    genome_length += line_len - 2;
    sequence = (char **) realloc (sequence, (genome_length + N_contigs) * sizeof(char));
    if ( sequence[N_contigs - 1] == NULL) {
        sequence[N_contigs - 1] = (char *) malloc (line_len * sizeof(char));
        strcpy(sequence[N_contigs - 1], line.c_str());
    }
    else {
        sequence[N_contigs -1] = (char *) realloc (sequence[N_contigs - 1], (strlen(sequence[N_contigs - 1]) + line_len) * sizeof(char));
        strcat(sequence[N_contigs - 1], line.c_str());
    }
    return 0;
}

int FastaParser::parse_fasta() {
    string line;
    int header_Nchar, line_len;
    header_Nchar = 0;
    if (fasta_file->good()) {
        getline( *fasta_file, line);
    }
    else {
        fprintf(stderr, "The fasta file cannot be parsed!\n");
        exit(0);
    }
    while ( fasta_file->good() ) {
        line_len = line.length() + 2;
        if (line.size() == 0) 
            ;
        else {
            if (line.at(0) == '>') {
                record_header( line, line_len, &header_Nchar);
                N_contigs++;
            }
            else {
                record_sequence( line, line_len );
            }
        }
        getline( *fasta_file, line );
    }
    return 0; 
}
