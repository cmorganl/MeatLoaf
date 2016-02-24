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
    int num_partitions;
    int total_seqs;
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

        num_partitions = 1;
        total_seqs = 0;
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
    \tThe output file. If not specified, output sequences in FASTA format are printed to stdout [OPTIONAL]\n\
    -m \n\
    \tThe minimum sequence length to be parsed from input [DEFAULT = 1]\n\
    -n: \n\
    \tNumber of partitions to divide the FASTA file.\n\
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
    while ((c = getopt (argc, argv, "i:o:n:m:hv")) != -1) {
        switch (c) {
            case 'i':
                this->input = optarg;
                break;
            case 'o':
                this->output = optarg;
                break;
            case 'n':
                this->num_partitions = atoi(optarg);
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

void write_subset_to_output(Fasta* genome, string output, int partition) {
    int seqs_written;
    int seqs_parsed = 0;
    int num_seqs_per = genome->N_contigs / partition;
    string header;
    std::ofstream fa_out;
    if (partition > genome->N_contigs) {
        std::cerr << "ERROR: Provided number of partitions is greater than number of sequences in input FASTA!" << std::endl;
        exit(-1);
    }

    while (partition >= 1) {
        seqs_written = 0;
        char str_partition [5];
        sprintf(str_partition, "%d", partition);
        string output_file = output + "_" + str_partition;
        fa_out.open(output_file.c_str());
        while (seqs_written <= num_seqs_per and seqs_parsed < genome->N_contigs) {
            header = genome->header_base.name[seqs_parsed].append("\n");
            fa_out.write(header.c_str(), header.length());
            fa_out.write(genome->sequences[seqs_parsed], genome->find_sequence_length(seqs_parsed));
            fa_out.write("\n", 1);
            seqs_written += 1;
            seqs_parsed += 1;
        }
        partition -= 1;
        fa_out.close();
    }
    return;
}

int main(int argc, char *argv[]) {
    Options options;
    options.Options::FetchArguments(argc, argv);
    if (options.hFlag == true) {
        options.Print_help();
        exit(1);
    }
    Fasta genome(options.input);

    int return_status = genome.parse_fasta(options.min_length);
    if (return_status > 0)
        fprintf(stderr, "ERROR: The input was not parsed correctly (N_contigs differs from vector size)!");

    std::cout << "Longest sequence = " << genome.find_longest_contig() << std::endl;
    std::cout << "Number of bases in file = " << genome.genome_length << std::endl;


    write_subset_to_output(&genome, options.output, options.num_partitions);


}