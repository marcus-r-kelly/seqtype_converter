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

    char *infilename = NULL;
    char *outfilename= NULL ; 
    char forcedIFF[40] ; 
    //char *inputFile ; 
    seqFormat inputFileFormat ; 
    seqFormat outputFileFormat ; 
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


    while ( (c = getopt (argc, argv, "Di:o:f:I:hH") ) != -1)
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
                infilename= new char[40] ; 
                strcpy(infilename,optarg) ; 
                break ;  
            case 'o':
                // determine input file name
                outfilename= new char[40] ; 
                strcpy(outfilename,optarg) ; 
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


    if ( argc > optind && infilename == NULL ) 
    {
        // if you haven't provided the filename and some string is still out there, assume that
        // the last argument is the input file name
        if ( debug ) cerr << "DEBUG: Looking for implicit input file name. (argc: " << argc << ")" << endl ; 
        infilename= new char[40] ; 
        strcpy(infilename,argv[optind]) ; 
    }
    else if ( infilename == NULL ) 
    {
        infilename= new char[40] ; 
        strcpy(infilename,"infile.txt") ; // a la phylip
    }

    if ( argc > optind && outfilename == NULL ) 
    {
        outfilename= new char[40] ; 
        strcpy(outfilename,"outfile.txt") ; // a la phylip
    }


    infile.open(infilename) ; 
    if (debug ) cerr << "File " << infilename << " opened for reading." << endl ; 
    

    if ( ! forceIFF)
    {
        inputFileFormat=getFileFormat(infile) ; 
        if ( debug ) cerr << "DEBUG: Interpreted format: " << printFileFormat(inputFileFormat) << endl ;
    }
    else
    {
        inputFileFormat=readFileFormat(forcedIFF)  ; 
        if (debug ) cerr << "Input format forced " << forceIFF << " as " << forcedIFF << endl ; 
    }

    switch (inputFileFormat)
    {
        case FASTA:
            theAln=readFASTA(infile ) ;
            break ; 
        case PHYLIP:
            theAln=readPHYLIP(infile ) ; 
            break ; 
        case PHYML:
            theAln=readPHYML(infile) ; 
            break ; 
        case GENBANK:
            theAln=readGENBANK(infile) ; 
            break ; 
        case CRAP:
        default:
            perror("Interpreted some kind of crappy file. Get it together.") ; 
            exit(EXIT_FAILURE);
    }

    if (debug ) cerr << "Output format: " << outputFormat << endl
    << "  Output file: " << outfilename << "opened for writing." << endl ; 
    outputFileFormat=readFileFormat(outputFormat) ; 

    outfile.open(outfilename) ; 

    switch ( outputFileFormat )
    {
        case FASTA:
            writeFASTA(outfile, theAln) ;
            break; 
        case PHYLIP:
            writePHYLIP(outfile, theAln) ;
            break; 
        case PHYML:
            perror("Not yet ready.") ; 
            exit(EXIT_FAILURE) ; 
            break ; 
        case GENBANK:
            perror("Not yet ready.") ; 
            exit(EXIT_FAILURE) ; 
            break ; 
        case CRAP:
        default:
            perror("Interpreted some kind of crappy file. Get it together.") ; 
            exit(EXIT_FAILURE);
    }

        



    return 0 ; 


}
