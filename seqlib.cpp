#include<string>
#include<iostream>
#include "seqlib.h" 
#include<cstdio>
#include<errno.h>
#include<cstdlib>
using namespace std ; 

const int DEFAULT_ALIGNMENT_SIZE=100 ; 

seq::seq(){
    size=0 ; 
    seqtype=DNA ; 
    name="" ; 
    contents="" ; 
}   

//  seq::seq(ifstream & file)
//  {
//      set(file) ; 
//  }

//  seq::seq(string fname)
//  {
//      set(fname) ; 
//  }

//  seq::~seq(){}  

void seq::setName(char* newname)
{
    string strname ;
    strname.append(newname) ; 
    this->name=strname; 
}

void seq::setName( string newname)
{
    this->name=newname;  
}

void seq::setContents( string str )
{
    if ( str.find(' ') != std::string::npos)// preliminary screen for spaces
    {
        deSpace(str) ; 
    }

    contents=str ; 
    size=str.length() ; 

    seqtype=determineType(str); 
}

void seq::setContents(char* str)
{
    contents.clear() ; 
    contents.append( str) ; 

    if ( this->contents.find(' ') != std::string::npos)// preliminary screen for spaces
    {
        deSpace(this->contents) ; 
    }

    size=contents.length() ; 

    seqtype=determineType(contents) ;
}


void seq::append(string str)
{

    // all kinds of things can go wrong here if str is empty
    if ( str.empty() )
    {
        cerr << "DEBUG: attempt to append empty string." << endl ; 
        return; 
    }

    string alphabet ;

    if ( str.find(' ') != std::string::npos)// preliminary screen for spaces
    {
        deSpace(str) ; 
    }

    if ( contents.empty() )
    {
        cerr << "DEBUG : Initializing type ; input appended to empty contents." << endl ; //DEBUG
        seqtype=determineType(str) ; 
        cerr << "DEBUG : determined type " << typeStr(seqtype) << endl ; 
    }

    switch ( this->seqtype ){
        case DNA:
            alphabet=DNA_ALPHABET ; 
        case RNA:
            alphabet=RNA_ALPHABET ; 
        case AA:
            alphabet=AA_ALPHABET ; 
        default:
            cerr << "Unrecognized sequence type, function append(string)" << endl ; 
            exit(EXIT_FAILURE) ; 
    }

    if ( contents != "" && ! testType(str,alphabet) )
    {

        cerr << "Input " << str << "produces mismatch between stored type " << this->getTypeStr() <<
        "and inferred type " << typeStr(determineType(str)) << endl ; 

        perror("Appended string does not match type of currently stored string") ;
        exit(1) ; 
    }

    contents.append(str) ; 
}

void seq::append(char* str)
{
    string cppstr ;
    cppstr.append(str) ; 
    string alphabet ; 

    if ( cppstr.empty() )
    {
        cerr << "DEBUG: attempt to append empty string." << endl ; 
        return; 
    }

    if ( cppstr.find(' ') != std::string::npos)// preliminary screen for spaces
    {
        deSpace(cppstr) ; 
    }


    if ( contents.empty() )
    {
        cerr << "DEBUG : Initializing type ; input appended to empty contents." << endl ; //DEBUG
        seqtype=determineType(cppstr) ; 
        cerr << "DEBUG : determined type " << typeStr(seqtype) << endl ; 
    }

    switch ( this->seqtype ){
        case DNA:
            alphabet=DNA_ALPHABET ;
            break ; 
        case RNA:
            alphabet=RNA_ALPHABET ; 
            break ; 
        case AA:
            alphabet=AA_ALPHABET ; 
            break ; 
        default:
            cerr << "Unrecognized sequence type, function append(char*)" << endl ; 
            exit(EXIT_FAILURE) ; 
    }

    if ( contents != "" && ! testType(cppstr,alphabet) )
    {
        cerr << "Input " << cppstr << " produces mismatch between stored type " << this->getTypeStr() <<
        " and inferred type " << typeStr(determineType(str)) << endl ; 
        perror("Appended string does not match type of currently stored string") ;
        exit(1) ; 
    }

    contents.append(cppstr) ; 
}

SeqType seq::determineType(string query)
{

    if ( testType(query,DNA_ALPHABET) )
    {
      //  seqtype=DNA ; 
            return DNA; 
    }
    else if ( testType(query,RNA_ALPHABET) ) 
    {
       // seqtype=RNA ; 
         return RNA; 
    }
    else if ( testType(query,AA_ALPHABET) ) 
    {
        //seqtype=AA ; 
        return AA ; 
    }
    else
    {
        cerr << "Unrecognizable sequence type." << endl ; 
        //seqtype=DNA; 
        return DNA; 
    }

}

bool seq::testType(std::string query,std::string alphabet)
{
    unsigned int i ; 

    for ( i=0; i < query.size() ; i++ ) 
    {
        if ( alphabet.find(query[i]) == string::npos ) 
        {
            cerr << "DEBUG: Offending character: " << query[i] << endl ; 
            return false ; 
        }
    }

    return true ; 
}

string seq::getName()
{
    return name ; 
}

string seq::getSeq()
{
    return contents ; 
}

int seq::length() 
{
    return size ; 
}

char & seq::operator[](int index)
{
    return contents[index] ; 
}

string seq::getTypeStr()
{
    switch(seqtype){
        case DNA:
            return "DNA" ;
            break ; 
        case RNA:
            return "RNA" ; 
            break ; 
        case AA : 
            return "AA" ; 
            break ; 
        default:
            cerr << "ERROR: getTypeStr(): invalid sequence type" << endl ; 
            exit(1); 
    }
}

SeqType seq::getType()
{
    return seqtype ; 
}



/*
 * END OF THINGS THAT RELATE TO SEQ, AND BEGINNING OF THINGS THAT RELATE TO ALN
 */

aln::aln()
{
    numTaxa=0 ; 
    alloc_taxa=16 ; 
    numChars=0 ; 
    alntype=DNA ; 

    names=new string[alloc_taxa] ; 
    sequences=new seq*[alloc_taxa] ; 
    cols=new int[alloc_taxa] ; 
}

aln::~aln()
{
    delete[] names ; 
    delete[] sequences ; 
    delete[] cols ; 
}

string aln::getNames()
{
    int i ;
    string out ; 

    for ( i=0 ; i < numTaxa ; i++ )
    {
        out += sequences[i]->getName() ; 
        out += " " ; 
   }

    out.erase(out.end()) ; 

    return out ; 
}

void aln::add(seq* newseq)
{
    if ( numTaxa == alloc_taxa) 
        resize() ; 

    if ( numTaxa==0 && alntype != newseq->getType() ) 
    {
        // the addition of the first sequence sets
        // the type for the whole alignment
        alntype=newseq->getType();
    }
    else if ( alntype != newseq->getType() ) 
    {
        cerr << "ERROR: aln::add() Clashing types: sequence " << newseq->getTypeStr()
             << " alignment: " << this->getTypeStr() << endl ; 
             exit(1) ; 
    }

    numTaxa += 1 ; 
    numChars+=newseq->length() ; 

    sequences[numTaxa-1]=newseq ; 
    names[numTaxa-1]=newseq->getName() ; 
    cols[numTaxa-1]=newseq->length() ; 

}

void aln::resize()
{
    int i ;
    string* newnames ; 
    seq** newseqs ; 
    int* newcols ;

    newnames=new string[2*alloc_taxa] ; 
    newseqs=new seq*[2*alloc_taxa] ; 
    newcols=new int[2*alloc_taxa] ; 

    for ( i=0 ; i < numTaxa ; i++ ) 
    {
        newnames[i]=names[i] ; 
        newseqs[i]=sequences[i] ; 
        newcols[i]=cols[i] ; 
    }

    alloc_taxa*=2 ; 

    delete[] names ; 
    delete[] sequences ; 
    delete[] cols ; 

    names=newnames ; 
    sequences=newseqs ; 
    cols=newcols ; 
}

string aln::getTypeStr()
{
    switch(alntype){
        case DNA:
            return "DNA" ;
            break ; 
        case RNA:
            return "RNA" ; 
            break ; 
        case AA : 
            return "AA" ; 
            break ; 
        default:
            cerr << "ERROR: getTypeStr(): invalid sequence type" << endl ; 
            exit(1); 
    }
}

SeqType aln::getType()
{
    return alntype ; 
}

int* aln::lengths()
{
    return cols ; 
}

bool aln::constantLength()
{
    int length ; 
    int i ;

    if ( numTaxa == 0 )
        return true ; 
    else
        length=sequences[0]->length() ; 


    for ( i=1 ; i < numTaxa ; i++ ) 
    {
        if ( sequences[i]->length() != length ) 
            return false ; 
    }

    return true ; 

}

int aln::taxa()
{
    return numTaxa ; 
}

int aln::chars()
{
    return numChars ; 
}

seq* & aln::operator[](int index)
{
    if ( index >= numTaxa ) 
    {
        perror("Index is greater than number of taxa.") ; 
        exit(1) ;
    }
    return sequences[index] ; 
}

seq* & aln::operator[](string name)
{
    int i=0 ; 
    for ( i=0 ; i < numTaxa ; i++ ) 
    {
        if ( names[i] == name )
        {
            return sequences[i] ; 
        }
    }

    cerr << "String " << name << "does not correspond to a valid sequence." << endl ; 
    perror("Invalid string index.") ; 
    exit (1) ; 
}

int aln::longest()
{
    int i ; 
    int soFar=0 ; 

    for ( i =0 ; i < this->taxa() ; i++ )
    {
        if ( cols[i] > soFar )
            soFar=cols[i] ;
    }

    return soFar ; 

}

int aln::shortest()
{
    int i ; 
    int soFar=-1 ; 

    for ( i =0 ; i < this->taxa() ; i++ )
    {
        if ( cols[i] < soFar || soFar==-1 )
            soFar=cols[i] ;
    }

    return soFar ; 

}

bool aln::uniform()
{
    int i ; 
    int soHi,soLo=-1 ; 

    for ( i =0 ; i < this->taxa() ; i++ )
    {
        if ( cols[i] < soLo || soLo==-1 )
            soLo=cols[i] ;
        if ( cols[i] > soHi )
            soHi=cols[i] ;
    }

    return ( soHi == soLo ) ;
   
}

// external functions provided by this library

string typeStr( SeqType seq_type)
{
    switch(seq_type){
        case DNA:
            return "DNA" ;
            break ; 
        case RNA:
            return "RNA" ; 
            break ; 
        case AA : 
            return "AA" ; 
            break ; 
        default:
            cerr << "ERROR: typeStr(): invalid sequence type" << endl ; 
            exit(1); 
    }
}


void deSpace(string & str)
{
    int i ; 
    for ( i=0 ; i <= (int) str.length() ; i++ )
    {
        if ( str[i] == ' ' ) // if you find a space
        {
            str.erase(str.begin()+i) ; // delete it
        }
    }

}
