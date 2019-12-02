#include "Fasta.hpp"

//#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>

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
    \tTurns on verbosity (prints summary ).\n\
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
    long seqs_written;
    long seqs_parsed = 0;
    long num_seqs_per = genome->N_contigs / partition;
    string header;
    string output_file;
    std::ofstream fa_out;
    if (partition > genome->N_contigs) {
        std::cerr << "ERROR: Provided number of partitions is greater"
                " than number of sequences in input FASTA!" << std::endl;
        exit(-1);
    }

    while (partition >= 1) {
        seqs_written = 0;
        if (partition >= 2) {
            char str_partition [5];
            sprintf(str_partition, "%d", partition);
            output_file = output + "_" + str_partition;
        }
        else
            output_file = output;

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

std::vector<long> sort_contigs_by_length(Fasta* seqs, bool verbose) {
    std::vector< long > seq_lengths (seqs->N_contigs);
    if (verbose)
        std::cout << "Sorting the sequences by length... " << std::flush;
    for (int x = 0; x < seqs->N_contigs; x++) {
        if ( seqs->sequences[x] == NULL) {
            fprintf(stderr, "ERROR: A sequence was improperly loaded into FASTA!\n");
            cout << "Index " << x << " (" << seqs->header_base.name[x] << ") is NULL." << endl;
            exit(-1);
        }
        seq_lengths.at(x) = strlen(seqs->sequences[x]);
        if (seq_lengths.at(x) != seqs->header_base.seq_length[x] )
            fprintf(stderr, "ERROR: Sequences and header_base are not matching!\n");
    }
    std::sort (seq_lengths.begin(), seq_lengths.end());

    if (verbose)
        std::cout << "done." << endl;
    return seq_lengths;
}


int median_seq_length(std::vector<long int> vector_lengths, long int num_seqs) {
    int half = num_seqs/2;
    return vector_lengths[half];
}


float calculateSD(std::vector<long int> data) {
    float sum = 0.0;
    float mean;
    float standardDeviation = 0.0;

    for (std::vector<long int>::iterator it = data.begin(); it != data.end(); ++it)
        sum += *it;
    mean = sum/data.size();

    for (std::vector<long int>::iterator it = data.begin(); it != data.end(); ++it)
        standardDeviation += pow(*it - mean, 2);

    return sqrt(standardDeviation / data.size());
}


int main(int argc, char *argv[]) {
    Options options;
    options.Options::FetchArguments(argc, argv);
    if (options.hFlag) {
        options.Print_help();
        exit(1);
    }
    Fasta genome(options.input);

    if (options.verbose)
        std::cout << "Reading " << options.input << "... " << std::flush;
    int return_status = genome.parse_fasta(options.min_length);
    if (options.verbose)
        std::cout << "done." << std::endl;

    if (return_status > 0)
        fprintf(stderr, "ERROR: The input was not parsed correctly (N_contigs differs from vector size)!");

    if (options.verbose) {
        std::cout << "Number of sequences = " << genome.N_contigs << endl;
        std::cout << "Number of bases in file = " << genome.genome_length << std::endl;
        std::vector<long int> sorted_seq_lengths = sort_contigs_by_length(&genome, options.verbose);
        float sd = calculateSD(sorted_seq_lengths);
        std::cout << "Longest sequence = " << sorted_seq_lengths.at( genome.N_contigs - 1) << std::endl;
        std::cout << "Shortest sequence = " << sorted_seq_lengths.at(0) << std::endl;
        std::cout << "Median length = " << median_seq_length(sorted_seq_lengths, genome.N_contigs) << std::endl;
        std::cout << "Sequence length standard deviation = " << sd << std::endl;
//        for (std::vector<int>::iterator it = sorted_seq_lengths.begin(); it != sorted_seq_lengths.end(); ++it)
//            std::cout << ' ' << *it;
//        std::cout << '\n';
    }


    write_subset_to_output(&genome, options.output, options.num_partitions);

}
