#include<iostream>
#include<fstream>
#include "seqlib.h"
#include<unistd.h>
#include<cctype>
#include<cstring>
#include "seqfilelib.h"
using namespace std ; 


bool debug=false ; 

void usage()
{
    cout << "ATM, you're SOL." << endl ; 
}


int main (int argc, char **argv)
{

    char *filename = NULL;
    char forcedIFF[40] ; 
    //char *inputFile ; 
    seqFormat inputFileFormat ; 
    bool forceIFF=false ; 

    ifstream infile ; 
    ofstream outfile ;
    int c;
    aln theAln ;

    char outputFormat[40] ; 

    strcpy(outputFormat,"phylip") ; 

    opterr = 0;

    if ( argc == 0)
    {
        usage() ; 
        return 0 ; 
    }


    while ( (c = getopt (argc, argv, "Di:f:I:hH") ) != -1)
    {
        switch (c)
        {
            case 'D':
                debug=true ;    
                if (debug ) cerr << "DEBUG: activated." << endl ; 
                break ; 
            case 'f':
                // determine output file format
                strcpy(outputFormat,optarg)  ;
                break;
            case 'I' :
                // forces input file format
                forceIFF=true ; 
                strcpy(forcedIFF,optarg)  ;
                break ; 
            case 'i':
                // determine input file name
                filename= new char[40] ; 
                strcpy(filename,optarg) ; 
                break ;  
            case 'H':
                // call for help
                usage() ;
                return 0 ; 
                break;
            case 'h':
                // call for help
                usage() ;
                return 0 ; 
                break;
            case '?':
                if (optopt == 'f')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                    return 1;
            default:
                break ; 
       }
    }


    if ( argc > optind && filename == NULL ) 
    {
        // if you haven't provided the filename and some string is still out there, assume that
        // the last argument is the input file name
        if ( debug ) cerr << "DEBUG: Looking for implicit input file name. (argc: " << argc << ")" << endl ; 
        filename= new char[40] ; 
        strcpy(filename,argv[optind]) ; 
    }
    else if ( filename == NULL ) 
    {
        filename= new char[40] ; 
        strcpy(filename,"infile.txt") ; // a la phylip
    }

    if (debug ) cerr << "Output format: " << outputFormat << "  Input file: " << filename << endl ; 

    if (debug && forceIFF ) cerr << "Input format forced " << forceIFF << " as " << forcedIFF << endl ; 

    infile.open(filename) ; 
    if (debug ) cerr << "File " << filename << " opened for reading." << endl ; 
    

    if ( ! forceIFF)
    {
        inputFileFormat=getFileFormat(infile) ; 
        if ( debug ) cerr << "DEBUG: Interpreted format: " << printFileFormat(inputFileFormat) << endl ;
    }

    switch (inputFileFormat)
    {
        case FASTA:
            theAln=readFASTA(infile ) ;
        case PHYLIP:
            theAln=readPHYLIP(infile ) ; 
        case PHYML:
            theAln=readPHYML(infile) ; 
        case GENBANK:
            theAln=readGENBANK(infile) ; 
        case CRAP:
        default:
            perror("Interpreted some kind of crappy file. Get it together.") ; 
            exit(1);
    }
        



    return 0 ; 


}
