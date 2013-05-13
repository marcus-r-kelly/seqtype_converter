#include<string>
#ifndef _seqlib_h
#define _seqlib_h

#include<iostream>
#include<cstdio>
#include<errno.h>
#include<cstdlib>
using namespace std ; 

enum SeqType { DNA, RNA, AA } ; 
// NB that all of these allow for IUPAC degeneracies
const string DNA_ALPHABET="atgcyrswkmbdhvnATGCYRSWKMBDHVN-" ; 
const string RNA_ALPHABET="augcyrswkmbdhvnAUGCYRSWKMBDHVN-" ; 
const string AA_ALPHABET="ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz-" ; 

// external function provided by this library
string typeStr( SeqType seq_type) ; 
void deSpace(string & str) ; 

class seq{

public: 

    seq() ; 
    // leave unmade 
 // THOUGHT OF THE DAY: LEAVE DATA ENTRY/PREPROCESSING TO EXTERNAL FUNCTIONS
 // seq(ifstream &file) ; 
 // // read from file
 // seq(string fname) ; 
 // // read from file named fname

    ~seq() ; 
    // destructor

    void setName(string name) ; // set name of sequence
    void setName(char* name) ;  // ^^
    void setContents(string str) ; // seq sequence contents to str
    void setContents(char* word) ; 
    void append(string str) ; // append str to sequence contents
    void append(char* str) ;  // ^^

    string getName() ;
    string getName(int start, int finish) ;
    // get name of single sequence

    string getSeq() ; 
    // get sequence
    string getSeq(int start, int finish) ;
    // get sequence data from indices start to finish, inclusive

    int length() ; // returns sequence length (var size)


    char& operator[](int index) ; // returns character at position index of sequence

    string getTypeStr() ; // returns sequence type as string
    SeqType getType() ; // return sequence type as SeqType

private:

    int size ; // stores sequence length

    SeqType seqtype ; //stores sequence type
    string name;     // stores sequence name
    string contents; // stores sequence contents

    SeqType determineType(string query) ; // initializes type from string contents
    bool testType(std::string query,std::string alphabet) ; // sees if types are valid

} ; 

class aln{

    public:

        aln() ; 
     // // leave usued
     // aln(ifstream &file) ; 
     // // from file
     // aln(string fname) ; 
     // // from file named fname

        ~aln() ; 

        string getNames() ; 
        //returns string of concatenated sequence names

        void add(seq* newseq) ; // adds newseq to sequence alignment

        int* lengths() ; // returns pointer to array of ints of string lenghts
        bool constantLength() ; // tests if all sequence strings are the same length

        int taxa() ;
        int chars() ; 

        seq* & operator[](int index) ; 
        seq* & operator[](string name) ; 

        string  getTypeStr() ; 
        SeqType getType() ; 

        int longest() ; 
        int shortest() ; 
        bool uniform() ; 

        bool uniqueNames() ; 
        bool uniqueNames(int nochars) ; 

    private:

        int* cols;

        int numTaxa ; 
        int alloc_taxa ;

        int numChars ;

        SeqType alntype ; 
        string* names;
        seq** sequences ; 

        void resize() ; 

} ;

#endif
