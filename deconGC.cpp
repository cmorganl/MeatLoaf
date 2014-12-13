#include <stdio.h>
#include <string.h>
#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <getopt.h>

using namespace std;

struct fasta {
    char *head;
    char *seq;
};

struct seq {
    char *header;
    int G;
    int C;
    double Total;
};

int checkArgs(int hflag, int argc, char *argv[]) {
    string USAGE = "USAGE: \t %s [options] -i assembly.fasta\n\n\
-h\t\tShow help message and exit\n\
-o\t\toutputPrefix[OPTIONAL]\n\
-d [DEFAULT] \tOR\t-v \n\n";
    int USAGE_len = USAGE.length();
    char *usage = (char*) malloc (USAGE_len * sizeof(char));
    strcpy(usage, USAGE.c_str());

    if (argc < 2) {
        printf("An assembly file is required, one was not provided.\n\n");
        printf(usage, argv[0]);
        exit(0);
    }
    else if (hflag == 1) {
        string HELP = "Help description:\n\
-i: \n \
\t The input file for GC analysis. This must be in FASTA format \n\
-o: \n \
\t The output prefix. Depending on which flag is selected \n\
\t (-d/-v this output file extension will change. \n\
-d: \n \
\t This will produce a fasta file in which the contigs with a GC content \n\
\t diverging by greater than 10\% of the mean GC content of the genome \n\
\t have been removed from the inputted FASTA file (decontaminated). \n\
\t This is the default operation. \n\
-v: \n \
\t A k-mer count file for all monomers (A, C, G, T) for every contig. \n";
       printf(usage, argv[0]);
       int HELP_len = HELP.length();
       char *help = (char*) malloc (HELP_len * sizeof(char));
       strcpy(help, HELP.c_str());
       printf("%s", help);
       exit(0);
    }
}

int decontaminate(struct fasta *genome, struct seq *seqs, int seqsN, int dflag, int vflag, char *output = NULL) {
    struct fasta *Filtered = (struct fasta*) malloc (seqsN * (Filtered, sizeof(fasta)));
    int n, i, header_len, contig_len, filteredN;
    header_len, contig_len, filteredN = 0;
    double GC_content, mean, grandGC = 0.0;
    for (n = 0; n < seqsN; n++) {
        GC_content = (((seqs[n].G + seqs[n].C)*100)/seqs[n].Total);
        grandGC += GC_content;
        if (vflag == 1)
            printf("%s : %.6f\n", seqs[n].header, GC_content);
        }
    mean = (grandGC/seqsN);

    printf("The mean GC content of the sample = %.2f\n", mean);
    printf("Flagged as possible contaminants:\n");
    for (i = 0; i<seqsN; i++) {
        GC_content = (((seqs[i].G + seqs[i].C)*100)/seqs[i].Total);
        header_len = strlen(seqs[i].header)+1;
        Filtered[i].head = (char*) malloc (header_len * sizeof(char));
        strcpy(Filtered[i].head, seqs[i].header);
        if (((GC_content + 10) < mean) || ((GC_content - 10) > mean)) {
            fprintf(stdout, "%s\t%.2f \n", seqs[i].header, GC_content);
            Filtered[i].seq = NULL;
        }
        else {
            contig_len = strlen(genome[i].seq)+1;
            Filtered[i].seq = (char*) malloc (contig_len * sizeof(char));
            strcpy(Filtered[i].seq, genome[i].seq);
        }
        filteredN++;

    }
    if (dflag == 1) {
        if (output != NULL) {
            FILE *decom_out;
            decom_out=fopen(output, "w");
            for (n = 0; n < filteredN; n++) {
                if (Filtered[n].seq != NULL) {
                    fprintf(decom_out, "%s\n", Filtered[n].head);
                    fprintf(decom_out, "%s\n", Filtered[n].seq);
                }
                else
                    ;
            }
            fclose(decom_out);
        }
        else
            fprintf(stderr, "An output file name is needed for the decontaminated genome to be written to.\n"); 
    }
    else {
        fprintf(stdout, "Completed successfully\n");
        exit(1);
    }
}

int getNucCounts(struct fasta *genome, int seqsN, int dflag, int vflag, char *output) {
    struct seq *seqs = (struct seq*) malloc (1*(seqs, sizeof(seq)));
    int x, i, contig_len, header_len;
    int G_count, C_count, total = 0;
    for (i = 0; i < seqsN; i++) {
        seqs = (struct seq*) realloc (seqs, (seqsN+1) * sizeof(struct seq));
        header_len = strlen(genome[i].head)+1;
        seqs[i].header = (char*) malloc (header_len * sizeof(char));
        strcpy(seqs[i].header, genome[i].head);

        contig_len = strlen(genome[i].seq);
        total = G_count = C_count = 0;
        for (x = 0; x < contig_len; x++) {
            if (genome[i].seq[x] == 'G') {
                G_count++;
                total++;
            }
            if (genome[i].seq[x] == 'C') {
                C_count++;
                total++;
            }
            else if (genome[i].seq[x] == 'A' || genome[i].seq[x] == 'T')
                total++;
            else
                ;
        }
        seqs[i].G = G_count;
        seqs[i].C = C_count;
        seqs[i].Total = total;
    }
/*
    for (x = 0; x < seqsN; x++) 
        printf("%s : %d %d %.2f\n", seqs[x].header, seqs[x].G, seqs[x].C, seqs[x].Total);
*/    
    decontaminate(genome, seqs, seqsN, dflag, vflag, output);
}

main(int argc, char *argv[]) {

    char *input = NULL;
    char *output = NULL;
    int output_len;
    int dflag = 0;
    int vflag = 0;
    int hflag = 0;
    int index;
    int c;
    opterr = 0;
    while((c = getopt (argc, argv, "i:o:hvd")) != -1) {
        switch (c) {
            case 'i':
                input = optarg;
                break;
            case 'h':
                hflag = 1;
                break;
            case 'd':
                dflag = 1;
                break;
            case 'v':
                vflag = 1;
                break;
            case 'o':
                output = optarg;
                break;
            case '?':
                if ((optopt == 'i') || (optopt == 'o')) {
                    fprintf(stderr, "Option -%c requires an argument. \n", optopt);
                    exit(1);
                }
                else if (isprint (optopt))
                    fprintf(stderr, "Unknown option character -%c. \n", optopt);
                else
                    fprintf(stderr, "Unknown option `\\x%x'.\n", optopt);
                return 1;
            default:
                exit(0);
        }
    }

    checkArgs(hflag, argc, argv);

    ifstream fastaFile (input);
    if ( !fastaFile.is_open()) {
        printf("Could not open file %s\n", input);
        exit(0);
    }
    struct fasta *genome = (struct fasta*) malloc (1*(genome, sizeof(fasta)));
    int x, line_len, seqsN, first;
    char *line = (char*) malloc (1 * sizeof(char));
    string s;
    x = line_len = seqsN = first = 0;
    while( fastaFile.good() ) {
        getline( fastaFile, s);
        line_len = s.length();
        line = (char*) realloc (line, (line_len+1) * sizeof(char));
        strcpy(line, s.c_str());
        if (line[0] == '>') {
            genome = (struct fasta*) realloc (genome, (seqsN+1) * sizeof(struct fasta));
            genome[seqsN].head = (char*) malloc ((line_len+1) * sizeof(char));
            strcpy(genome[seqsN].head, line);
            first = 0;
            seqsN++;
        }
        else {
            if (first == 0) {
                genome[seqsN-1].seq = (char*) malloc ((line_len+1) * sizeof(char));
                strcpy(genome[seqsN-1].seq, line);
                first = 1;
            }
            else {
                genome[seqsN-1].seq = (char*) realloc (genome[seqsN-1].seq, (line_len + strlen(genome[seqsN-1].seq)+1) * sizeof(char));
                strcat(genome[seqsN-1].seq, line);
            }
        }
    }
    fastaFile.close();
    free(line);
    getNucCounts(genome, seqsN, dflag, vflag, output);
}
