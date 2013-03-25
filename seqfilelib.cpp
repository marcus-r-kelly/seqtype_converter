#include<string>
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<errno.h>
#include<fstream>
#include "seqfilelib.h"
using namespace std;

seqFormat getFileFormat(ifstream & infile )
{
    int i ; 
    int j ;
    int lastNonAlpha=0 ; 
    int maxLNA=0 ; 
    char line[MAX_LINE_SIZE] ; 
    char* word ; 
    int  phylip_taxa,phylip_length=-1 ; 

    infile.getline(line,MAX_LINE_SIZE) ;

    // First, we attempt to process any file as though
    // it contains the PHYLIP headers of taxa number and
    // alignment length.
    word=strtok(line," \t") ; 
    phylip_taxa=atoi(word) ;
    word=strtok(NULL," \t") ; 
    phylip_length=atoi(word) ; 

    if ( phylip_taxa > 0  && phylip_length > 0 )
    {
         // having determined that the sequence is of a phylip-ish nature,
        // we must down determine whether it adheres to more formal PHYLIP rules.
        // specifically, how many characters in is the last non-alphanumeric character
        // of ANY line?
        // The idea is that some character from some taxa name will occur after the 10th
        // character, which is allowed under PhyML format but NOT allowed under PHYLIP format.

        // Finally, it's worth noting that this is NOT a foolproof test, and in fact the
        // exact reason why users are allowed to force input file formats.

        i=0 ; 
        while ( infile.good() && i < phylip_taxa )
        {
            infile.getline(line,MAX_LINE_SIZE) ; 
            for ( j=0 ;  line[j] != '\0' ;  j++ ) 
            {
                if ( ! isalpha(line[j]) && line[j] != '-' && line[j] != ' ')
                // dashes in alignment gaps shouldn't be counted for this purpose
                // for readability purposes, phylip format will sometimes put spaces in
                {
                    lastNonAlpha=j ;
                }
            }

                cerr << "Last NA on line " << i << ":" << line[j] << endl ; 

            if ( lastNonAlpha > maxLNA)
            {
                maxLNA=lastNonAlpha ; 
            }

            i++ ;
        }

        if ( maxLNA < 10 )
        {
            cerr << "maxLNA: " << maxLNA << endl ; 
            return PHYLIP ; 
        }
        else if ( maxLNA < 20 )
        {
            cerr << "maxLNA: " << maxLNA << endl ; 
            return PHYML ; 
        }
        else
        {
            cerr << "Error finding last non-alphanumeric character in PHYLIP-like alignment"
                 << " body. (position " << maxLNA << "). Defaulting to PhyML." << endl;
                return PHYML ; 
        }
    }

    // otherwise, we try to differentiate between FASTA and genbank

    infile.seekg(ios::beg) ; 
    while ( infile.good() ) 
    {
        infile.getline(line,MAX_LINE_SIZE) ; 

        // signals a name of a sequence. If it starts a line, this is a giveaway for FASTA
        if ( strncmp(line,">",1) == 0 )
        {
            return FASTA ; 
        }
        // genbank formats are littered with these pdb-like headers
        else if ( strncmp(line,"FEATURES",8) == 0 )
        {
            return GENBANK ; 
        }
        else if ( strncmp(line,"COMMENT",7) == 0 )
        {
            return GENBANK ; 
        }
        else if ( strncmp(line,"BASE COUNT",10) == 0 )
        {
            return GENBANK ; 
        }
        else if ( strncmp(line,"ORIGIN",6) == 0 )
        {
            return GENBANK ; 
        }
        else if ( strncmp(line,"SOURCE",6) == 0 )
        {
            return GENBANK ; 
        }

    }

    perror("Stumped.") ; 
    return CRAP;

}

string    printFileFormat(seqFormat format) 
{
    string out ; 
    switch(format)
    {
        case FASTA:
            out.append("FASTA") ; 
            break ; 
        case PHYLIP:
            out.append("PHYLIP") ; 
            break ; 
        case PHYML:
            out.append("PHYML") ;
            break ; 
        case GENBANK:
            out.append("GENBANK") ; 
            break ;
        case CRAP:
            out.append("CRAP") ; 
            break ; 
    }

    return out ;

}

seqFormat readFileFormat(string format) 
{

    if ( format == "fasta" )
        return FASTA ; 
    if ( format == "FASTA" )
        return FASTA ; 
    if ( format == "phylip" )
        return PHYLIP ; 
    if ( format == "PHYLIP" )
        return PHYLIP ; 
    if ( format == "phyml" )
        return PHYML;
    if ( format == "PHYML" )
        return PHYML;
    if ( format == "PhyML" )
        return PHYML;
    if ( format == "genbank" )
        return GENBANK  ;
    if ( format == "GENBANK" )
        return GENBANK  ;
    if ( format == "GenBank" )
        return GENBANK  ;
    else
    {
        perror("String could not be interpreted as file format. Exiting") ; 
        exit(1) ; 
    }

}

seqFormat readFileFormat(char* format) 
{
    string formatstring ; 
    formatstring.append(format) ; 


    return readFileFormat(formatstring) ; 

}


aln & readFASTA(ifstream & infile ) 
{

    aln* newAln = new aln ; 
    char line[MAX_LINE_SIZE] ;
    seq* currentSeq=NULL ; 
    string currSeqName ; 

    while ( infile.good() )
    {
        infile.getline(line,MAX_LINE_SIZE) ; 

        if ( line[0] == '>' )
        {
            if ( currentSeq !=NULL ) 
            {
                newAln->add(currentSeq) ; 
            }
            currentSeq=new seq ; 
            currentSeq->setName(line + 1) ; 
            cerr << "Detected new sequence named " << currentSeq->getName() << endl ; 
        }
        else
        {
            currentSeq->append(line) ; 
        }
    }

    cerr << "Contained " << newAln->taxa() << " taxa." << endl ;  
    cerr << "Sequence lengths:" << endl ; 
    for ( int i=0 ; i < newAln->taxa() ; i++ )
    {
        cerr << (newAln->lengths())[i] << " "  ;
        if ( i % 5 == 0 )
            cerr << endl ; 
    }

    return *newAln ; 

}

aln & readPHYLIP(ifstream & infile)
{

    aln* newAln = new aln ; 
    char    line[MAX_LINE_SIZE] ;
    char    word[MAX_LINE_SIZE] ; 
    seq*    currentSeq=NULL ; 
    string  currSeqName ; 

    int phylip_taxa,phylip_length=1 ; 
    unsigned int blockLength=0;
    int linesReadThisBlock=0 ;
    int blocksRead=0 ;
    int totalRead=0 ; 
    int lineno=1 ; 
    bool firstBlock=true ; 

    infile.getline(line,MAX_LINE_SIZE) ;

    word=strtok(line," \t") ; 
    phylip_taxa=atoi(word) ;
    word=strtok(NULL," \t") ; 
    phylip_length=atoi(word) ; 

    // this loop begins on the second line of the file

    while ( infile.good() && totalRead < phylip_taxa*phylip_length ) 
    {
        infile.getline(line,MAX_LINE_SIZE) ; 

        if (strlen(line) < 10 ) // to counter for blank lines
            continue ; 

        if ( linesReadThisBlock == phylip_taxa - 1) 
        {
            firstBlock=false ;
            linesReadThisBlock=0 ; 
            blocksRead++ ; 
        }


        if ( firstBlock )
        {
            // copy over name
            strncpy(word,line,10) ; 
            word[10]='\0' ;
            currentSeq=new seq ;
            currentSeq->setName(word) ; // on the first block, set the sequence name

            cerr << "Initialized sequence " << linesReadThisBlock << " with name " << word << "." << endl ; 

            strncpy(word,line + 10,MAX_LINE_SIZE) ;
            /*if ( blockLength == 0 ) 
                blockLength=strlen(word) ;  // initialize block length if it hasn't been set already
            else if (blockLength != strlen(word) ) ; // something really needed to go wrong to set this off
            {
                cerr << "Inconsistent block lengths at line " << lineno << ":" << endl << line << endl ; 
                cerr << "Previous block length was " << blockLength << ", current is " << strlen(word) << endl ;  
                exit(1) ; 
            }*/

            currentSeq->setContents(word) ;
            // we dont' need to worry about removing arbitrary spaces because
            // the seq data type now checks for them each time its contents are updated.

            // at this point, currentSeq has been created and contains the first block's worth
            // of sequence data. it has NOT been assigned to an alignment.

            newAln->add(currentSeq) ; 

            cerr << "Appended " << strlen(word) << " characters (including spaces) to sequence " <<
                currentSeq->getName() << "." << endl ; 

            linesReadThisBlock++ ; 
            totalRead += strlen(word) ;

            cerr << "Total read: " << newAln.chars() << endl ; 

        }
        else
        {

            currentSeq=newAln[linesReadThisBlock] ;

            cerr << "Appended " << strlen(word) << " characters (including spaces) to sequence " <<
                currentSeq->getName() << "." << endl ; 

            strncpy(word,line+10,MAX_LINE_SIZE) ;
            currentSeq.append(word) ;

            linesReadThisBlock++ ; 
            totalRead += strlen(word) ;

            cerr << "Total read: " << newAln.chars() << endl ; 
        }
        lineno++ ; 
    }

    if ( newAln.chars() != phylip_taxa*phylip_length )
    {
        cerr << "ERROR: Read " << newAln.chars() << " characters rather than expected "
             << phylip_taxa*phylip_length << " characters." << endl ; 

        exit(1) ; 
    }


    return newAln ;

}

aln & readPHYML(ifstream & infile)
{

    aln* newAln = new aln ; 
    char    line[MAX_LINE_SIZE] ;
    char    word[MAX_LINE_SIZE] ; 
    seq*    currentSeq=NULL ; 
    string  currSeqName ; 

    int phyml_taxa,phyml_length=1 ; 
    unsigned int blockLength=0;
    int linesReadThisBlock=0 ;
    int blocksRead=0 ;
    int totalRead=0 ; 
    int lineno=1 ; 
    bool firstBlock=true ; 

    infile.getline(line,MAX_LINE_SIZE) ;

    word=strtok(line," \t") ; 
    phyml_taxa=atoi(word) ;
    word=strtok(NULL," \t") ; 
    phyml_length=atoi(word) ; 

    // this loop begins on the second line of the file
    // in what I'm calling "phyml" format, the 10-character name space has been relaxed.

    while ( infile.good() && totalRead < phyml_taxa*phyml_length ) 
    {
        infile.getline(line,MAX_LINE_SIZE) ; 

        if (strlen(line) < 10 ) // to counter for blank lines
            continue ; 

        if ( linesReadThisBlock == phyml_taxa - 1) 
        {
            firstBlock=false ;
            linesReadThisBlock=0 ; 
            blocksRead++ ; 
        }


        if ( firstBlock )
        {
            // copy over name
            word=strtok(line," \t") ; 

            currentSeq=new seq ;
            currentSeq->setName(word) ; // on the first block, set the sequence name

            cerr << "Initialized sequence " << linesReadThisBlock << " with name " << word << "." << endl ; 

            word=strtok(line," \n\t") ; 
            /*if ( blockLength == 0 ) 
                blockLength=strlen(word) ;  // initialize block length if it hasn't been set already
            else if (blockLength != strlen(word) ) ; // something really needed to go wrong to set this off
            {
                cerr << "Inconsistent block lengths at line " << lineno << ":" << endl << line << endl ; 
                cerr << "Previous block length was " << blockLength << ", current is " << strlen(word) << endl ;  
                exit(1) ; 
            }*/

            currentSeq->setContents(word) ;
            // we dont' need to worry about removing arbitrary spaces because
            // the seq data type now checks for them each time its contents are updated.

            // at this point, currentSeq has been created and contains the first block's worth
            // of sequence data. it has NOT been assigned to an alignment.

            newAln->add(currentSeq) ; 

            cerr << "Appended " << strlen(word) << " characters (including spaces) to sequence " <<
                currentSeq->getName() << "." << endl ; 

            linesReadThisBlock++ ; 
            totalRead += strlen(word) ;

            cerr << "Total read: " << newAln.chars() << endl ; 

        }
        else
        {

            currentSeq=newAln[linesReadThisBlock] ;

            cerr << "Appended " << strlen(word) << " characters (including spaces) to sequence " <<
                currentSeq->getName() << "." << endl ; 

            word=strtok(line," \n\t") ; 
            currentSeq.append(word) ;

            linesReadThisBlock++ ; 
            totalRead += strlen(word) ;

            cerr << "Total read: " << newAln.chars() << endl ; 
        }
        lineno++ ; 
    }

    if ( newAln.chars() != phyml_taxa*phyml_length )
    {
        cerr << "ERROR: Read " << newAln.chars() << " characters rather than expected "
             << phyml_taxa*phyml_length << " characters." << endl ; 

        exit(1) ; 
    }

    return newAln ;

}

seq& readGENBANK(ifstream & infile) ; 
{
    char    line[MAX_LINE_SIZE] ;
    char    word[MAX_LINE_SIZE] ; 
    seq*    currentSeq=new seq ; 

    infile.getline(line,MAX_LINE_SIZE) ;
    word=strtok(line," \t") ; 
    word=strtok(NULL," \t") ; 

    currentSeq.setName(word) ; // the accession number becomes the name of the sequence . 
    cerr << "Assigned sequence name " << currentSeq.getName() << endl ; 

    do
    {
        infile.getline(line,MAX_LINE_SIZE) ;
    } while ( strncmp(line,"ORIGIN",6) != 0 )

    do
    {
        word=strtok(line," \t") ; // hopefully iterates past the number markers
        currentSeq.append(word) ; 

    } while (strncmp(line,"//",2) != 0)
    // <-- marks termination of genbank files

    return currentSeq ; 

}

seq & readFASTA(ifstream & infile ) 
{

    char line[MAX_LINE_SIZE] ;
    seq* currentSeq=new seq ; 

    while ( infile.good() )
    {
        infile.getline(line,MAX_LINE_SIZE) ; 

        if ( line[0] == '>'  && (currentSeq->getName()).empty() )
        {
            currentSeq->setName(line+1) ; 
        }
        else if (line[0] == '>')
        {
            return *newAln ;  // read until second sequence
        }
        else
        {
            currentSeq->append(line) ; 
        }
    }

    return *newAln ; // or end of file

}

void writeFASTA( ofstream & outfile, aln & thealn)
{

    int seqno=0 ;
    int charno=0 ; 

    for ( seqno=0 ; seqno < thealn.taxa()  ; seqno++)
    {
        outfile << ">" << thealn[seqno]->getName() << endl ;
       
        charno=0 ; 
        while ( charno < thealn[seqno]->length() )
        {
            outfile.put( (thealn[seqno])[charno] ) ;

            if ( charno % FASTA_BLOCK_WIDTH == 0 )
                outfile.put('\n') ; 

        }

    }

} 

