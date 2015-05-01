/*
    This standalone RPKM calculator script was built upon 
    Evan Durno's countHits.cpp script from 2014. It has been modified
    to calculate RPKM for each sequence from the fasta file.
*/

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <getopt.h>

#include "FastaParser.hpp"
#include "Genome.hpp"

using namespace std ; 

struct multireads {
    char* read;
    char** sequences;
    int num_multi; //How many other sequences this read aligned to
};

struct alignmentsDict {
    char** sequences;
    float* num_alignments;
};
/*
bool isString( string a , string b )
{
	if( a.length() != b.length() )
		return false ; 
	for( int i = 0 ; i < a.length() ; i++ )
	{
		if( a[i] != b[i] )
			return false ; 
	}
	return true ; 
}

int findString( vector<string>* v , string str )
{
	for( int i = 0 ; i < v->size() ; i++ )
	{
		if( isString( v->at(i) , str ) )
			return i ; 
	}
	return -1 ; 
}
*/
int countHits( string samFile , char** headers, char** sequences, int N_contigs ) {
/*
Parsing the SAM file will occur in a single pass. Lines will be filtered for alignments,
indicated by a !* in the third field.
If there is an alignments, the bitflag is interrogated to determine whether the read was aligned to a 
single sequence or whether it was aligned to multiple sequences.
If the read was aligned once the sequence's accumulator will be increased in alignmentsDict.
If the read was aligned multiple times it will be added to multireads along with extra information.
*/
/*
	int N = orfs.size() ; 
	int out[N] ; 
	for( int i = 0 ; i < N ; i++ )
		out[i] = 0 ;
*/	
	string s;
	ifstream file( samFile.c_str() ) ; 
/*
Fields:
QNAME, bitFLAG, RNAME, POSition, MAPQuality
*/
	while( file.good() ) {
		getline( file , s ) ;
        while ( s[0] == '@' ) {
            getline( file, s );
        }
        if( s.length() > 0 ) { 
            cout << s << endl;
        }
//			if( s[0] != '@' ) {
//				int start ; 
//				for( start = 0 ; start < s.length() && s[start] != '\t' ; start++ ) { /* do nothing */	}
//				
//				start++ ; 
//				int end ; 
//				for( end = start ; end < s.length() && s[end] != '\t' ; end++ ) { /* do nothing */ } 
//				
//				end++ ; 
//				int val = atoi( s.substr( start , end - start ).c_str() ) ; 
//				
//				for( start = end - 1 ; start < s.length() && s[start] != '\t' ; start++ ) { /* do nothing */ } 
//				start++ ; 
//				for( end = start ; end < s.length() && s[end] != '\t' ; end++ ) { /* do nothing */ } 
//				end++ ; 
//				
//				string name = s.substr( start , end - start - 1 ) ; 
//				
//				if( ( val & ( 1 << four ) ) >> four  == 0 ) // if so, find name and increment 
//				{
//					int idx = findString( &orfs , name ) ; 
//					if( idx >= 0 )
//						out[idx]++ ; 
//				} 
//				n++ ; 
//			}
//		}
	}
	file.close() ; 
/*	
	for( int i = 0 ; i < N ; i++ )
		cout << '\t' << out[i] ; 
*/	
	return( 1 ) ; 
}

int main( int argc , char** argv ) {
/*
Inputs are a list of absolute paths to alignment files in SAM format
and a FASTA formatted file of sequences.

Output is a comma-separated matrix with columns being the sample name from the SAM file
and the FASTA headers are the rows.
*/
	if( argc < 3 ) { //Replace with getopt in the future
		printf("Requires a list of SAM files (argument 1) and a FASTA formatted file (argument 2)\n"); 
		printf("Example:\n\t./%s SAM_files.list contigs.fasta\n", argv[0]); 
		return( 1 ) ; 
	}

    FastaParser genome(argv[2]);
    genome.parse_fasta();
    cout << genome.N_contigs << endl;

	ifstream samList( argv[1] ); 
	string samFile; 
	while( samList.good() ) {
		getline( samList , samFile ) ;
		if( samFile.length() > 0 ) {
			cout << samFile ;
			countHits( samFile , genome.header, genome.sequence, genome.N_contigs ) ;
			cout << endl ;
		}
	}
	
	exit( 1 ) ; 
}
