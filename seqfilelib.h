
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

enum seqFormat {FASTA, PHYLIP, PHYML, GENBANK, CRAP} ; 

seqFormat getFileFormat(ifstream & infile ) ; 
string    printFileFormat(seqFormat format) ; 
seqFormat readFileFormat(string format) ; 


aln&    readFASTA(ifstream & infile ) ;
aln&    readPHYLIP(ifstream & infile ) ; 
aln&    readPHYML(ifstream & infile) ; 

seq&    readGENBANK(ifstream & infile) ; 
seq&    readFASTA(ifstream & infile ) ;
seq&    readPHYLIP(ifstream & infile ) ; 
seq&    readPHYML(ifstream & infile) ; 

aln&    writeFASTA(ifstream & infile ) ;
aln&    writePHYLIP(ifstream & infile ) ; 
aln&    writePHYML(ifstream & infile) ; 


#endif
