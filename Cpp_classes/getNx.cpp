#include "Fasta.hpp"

//#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <fstream>

struct Options {
// Fasta files:
    std::string input;
    std::string output;

// Parameters:
    int min_length;
    int errFlag;
    int opterr;
    bool problem;
    bool verbose;
    bool hFlag;

// Constructor:
    Options() {
        input = "";
        output = "";

        min_length = 1;

        errFlag = 0;
        opterr = 0;

        hFlag = false;
        verbose = false;
        problem = false;
    };
// Deconstructor:
    ~Options() {
        ;
    }
    void Print_usage(char *);
    void Print_help();
    int FetchArguments(int argc, char* argv[]);
};


void Options::Print_help() {
    std::string help;
    help = "\nHelp description:\n\
    -i: \n \
    \tThe FASTA input \n\
    -o \n\
    \tThe output file. If not specified, csv is printed to stdout [OPTIONAL]\n\
    -m \n\
    \tThe minimum sequence length to be parsed from input [DEFAULT = 1]\n\
    -v: \n\
    \tTurns on verbosity.\n\
    -h \n\
    \tShow help message and exit\n";
    std::cerr << help << std::endl;
    return;
}

void Options::Print_usage(char* arg) {
    std::cerr << "\nUSAGE: \t" << arg << " [options] -i input.fasta" << std::endl;
    return;
}

int Options::FetchArguments(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Type '" << argv[0] << " -h' for help options." << std::endl;
        this->problem = true;
        return 1;
    }
    int c;
    unsigned int j;
    char const * arguments_need_input = "i";
    while ((c = getopt (argc, argv, "i:o:m:hv")) != -1) {
        switch (c) {
            case 'i':
                this->input = optarg;
                break;
            case 'o':
                this->output = optarg;
                break;
            case 'm':
                this->min_length = atoi(optarg);
                break;
            case 'v':
                this->verbose = true;
                break;
            case 'h':
                this->hFlag = true;
                break;
            case '?':
                for (j=0; j < strlen(arguments_need_input); j++) {
                    if (optopt == arguments_need_input[j]) {
                        fprintf(stderr, "Option -%c requires an argument. \n", optopt);
                        exit(1);
                    }
                }
                if (isprint (optopt)) {
                    fprintf(stderr, "Unknown option character -%c. \n", optopt);
                    this->problem = true;
                }
                else {
                    fprintf(stderr, "Unknown option `\\x%x'.\n", optopt);
                    this->problem = true;
                }
                return 1;
            default:
                exit(0);
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    Options options;
    int return_status;
    return_status = options.Options::FetchArguments(argc, argv);
    if ( return_status > 0 )
        exit(1);
    if (options.hFlag == true) {
        options.Print_help();
        exit(2);
    }
    Fasta genome(options.input);

    return_status = genome.parse_fasta(options.min_length);
    if (return_status > 0)
        fprintf(stderr, "ERROR: The input was not parsed correctly (N_contigs differs from vector size)!");

    std::cout << "Number of sequences = " << genome.N_contigs << std::endl;
    std::cout << "Number of bases in file = " << genome.genome_length << std::endl;

    return_status = genome.writeNx(options.output, options.verbose);
    if (return_status != 1)
        exit(2);


}
