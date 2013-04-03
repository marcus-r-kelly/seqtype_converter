
#include<string>
#include<iostream>
#include<fstream>
#ifndef _seqfilelib_h
#define _seqfilelib_h
#include<cstdio>
#include<errno.h>
#include<cstdlib>
#include "seqlib.h"
using namespace std; 

const int MAX_LINE_SIZE=256 ; 
const int MAX_WORD_SIZE=16 ; 
const int FASTA_BLOCK_WIDTH=80 ; 
const int PHYLIP_BLOCK_WIDTH=60 ; 

enum seqFormat {FASTA, PHYLIP, PHYML, GENBANK, CRAP} ; 

seqFormat getFileFormat(ifstream & infile ) ; 
string    printFileFormat(seqFormat format) ; 
seqFormat readFileFormat(string format) ; 
seqFormat readFileFormat(char* format) ; 


aln&    readFASTA(ifstream & infile ) ;
aln&    readPHYLIP(ifstream & infile ) ; 
aln&    readPHYML(ifstream & infile) ; 

aln&    readGENBANK(ifstream & infile) ; 
//seq&    readFASTA_single(ifstream & infile ) ;

void   writeFASTA(ofstream & outfile, aln & thealn ) ;
void   writePHYLIP(ofstream & outfile, aln & thealn ) ; 
void   writePHYML(ifstream & outfile, aln & thealn) ; 


#endif
