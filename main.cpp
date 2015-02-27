#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include "version.h"

typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;


struct Options {
    int upper;
    int lower;
    int endLength;
    seqan::CharString bamFile;
    seqan::CharString baiIndexFile;
    seqan::CharString referenceFile;

    Options() : upper(-1), lower(3), endLength(500)
    {}
};
enum ContigEnd_t {ContigStart, ContigEnd};
typedef std::pair<std::pair<seqan::String<char>,ContigEnd_t>, std::pair<seqan::String<char>,ContigEnd_t> > MappingKey_t;
typedef std::map<MappingKey_t, int> Mapping_t;

class Graph {
    public:
        void addLink(seqan::String<char> id1, seqan::String<char> id2, ContigEnd_t e1, ContigEnd_t e2);
        void removeLinks(int lower = 3, int upper = -1);
    private:
        Mapping_t adjacency_list;

        friend std::ostream& operator<<(std::ostream& os, Graph& g);

};
void Graph::addLink(seqan::String<char> id1, seqan::String<char> id2, ContigEnd_t e1, ContigEnd_t e2) {
    MappingKey_t key;
    if( id1 < id2) {
        key = std::make_pair(std::make_pair(id1,e1), std::make_pair(id2,e2));
    } else {
        key = std::make_pair(std::make_pair(id2,e2), std::make_pair(id1, e1));
    }
    if(adjacency_list.find(key) == adjacency_list.end()) {
        adjacency_list[key] = 1;
    } else {
        adjacency_list[key] += 1;
    }
}

void Graph::removeLinks(int lower, int upper) {
    Mapping_t::iterator iter;
    for (iter = adjacency_list.begin(); iter != adjacency_list.end();) {
        if (lower >= 0) {
            if(iter->second < lower) {
                adjacency_list.erase(iter++);
                continue;
            }
        }

        if (upper >= 0) {
            if(iter->second > upper) {
                adjacency_list.erase(iter++);
                continue;
            }
        }
        ++iter;
    }
}

std::ostream& operator << (std::ostream& os, Graph& g) {
    Mapping_t::iterator iter;
    for(iter = g.adjacency_list.begin(); iter != g.adjacency_list.end(); ++iter) {

        os <<iter->first.first.first<<"\t"<<iter->first.second.first;
        os<<"\t"<<iter->second<<"\t";
        // first contig joined at begining or end?
        if(iter->first.first.second == ContigStart) {
            os << "start ";//"<-- ";
        } else {
            os << "end ";//"--> ";
        }

        if(iter->first.second.second == ContigStart) {
            os << "start";//"<--";
        } else {
            os << "end";//"-->";
        }
        os <<std::endl;
    }
    return os;
}

int parseRegion(int rID, int beginPos, int endPos, TBamContext const & context, seqan::BamIndex<seqan::Bai>& baiIndex, seqan::BamFileIn& inStream, Graph& g ) {
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inStream, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        //std::cerr << "ERROR: Could not jump to " << argv[3] << ":" << beginPos << "\n";
        return 1;
    }
    if (!hasAlignments)
        return 0;  // No alignments here.

    // Seek linearly to the selected position.
    seqan::BamAlignmentRecord record;
    while (!atEnd(inStream) )
    {
        readRecord(record, inStream );

        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || record.rID != rID || record.beginPos >= endPos)
            break;
        // If we are left of the selected position then we skip this record.
        if (record.beginPos < beginPos)
            continue;

        // make sure that the mate is not on the contig
        if(record.rID != record.rNextId && record.rNextId != -1) {
            // now check to see if we are mapping to an end of the other contig with that mate
            ContigEnd_t e1, e2;
            if(record.pNext < 500) {
                // maps to the begging of the contig
                e2 = ContigStart;
            } else if (record.pNext > seqan::contigLengths(context)[record.rNextId] - 500) {
                // maps to the end of the contig
                e2 = ContigEnd;
            } else {
                continue;
            }

            if(beginPos == 0) {
                // beggining of this contig
                e1 = ContigStart;
            } else {
                // end of this contig
                e1 = ContigEnd;
            }
            // add it to the graph.
            g.addLink(seqan::contigNames(context)[rID], seqan::contigNames(context)[record.rNextId], e1, e2);
        }
    }
    return 0;
}

seqan::ArgumentParser::ParseResult
parseCommandLine(Options& opts, int argc, char const ** argv) {
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bam2graph");
    seqan::setVersion(parser, PACKAGE_VERSION);
    seqan::setDate(parser, PACKAGE_DATE);

    seqan::addOption(parser, seqan::ArgParseOption( "l", "lower", "Lower coverage bound for the number of links between to contigs", seqan::ArgParseArgument::INTEGER, "INT")); 
    seqan::setDefaultValue(parser, "lower", "3");
    seqan::addOption(parser, seqan::ArgParseOption( "u", "upper", "Upper coverage bound for the number of links between to contigs", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "upper", "-1");
    seqan::addOption(parser, seqan::ArgParseOption( "e", "end-length", "The length from either side of the contig to search for links", seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "end-length", "500");
    seqan::addOption(parser, seqan::ArgParseOption( "b", "bam", "Bam file", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    seqan::setRequired(parser, "bam");
    seqan::addOption(parser, seqan::ArgParseOption( "B", "bai", "Index file", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    seqan::setRequired(parser, "bai");
    seqan::addOption(parser, seqan::ArgParseOption( "r", "ref-seqs", "Contigs to look for", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    seqan::setRequired(parser, "ref-seqs");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    //extract option values
    seqan::getOptionValue(opts.lower, parser, "lower");
    seqan::getOptionValue(opts.upper, parser, "upper");
    seqan::getOptionValue(opts.bamFile, parser, "bam");
    seqan::getOptionValue(opts.baiIndexFile, parser, "bai");
    seqan::getOptionValue(opts.referenceFile, parser, "ref-seqs");
    seqan::getOptionValue(opts.endLength, parser, "end-length");

    return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char const ** argv)
{
    
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Open Bam file for reading.
    seqan::BamFileIn input_bam(seqan::toCString(options.bamFile));
    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (!open(baiIndex, seqan::toCString(options.baiIndexFile)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << options.baiIndexFile << "\n";
        return 1;
    }

    seqan::BamHeader header;
    seqan::readHeader(header, input_bam);

    // get the contig names and lengths from the bam file
    TBamContext const & bamContext = context(input_bam);

    // read the references file
    std::ifstream refFile(seqan::toCString(options.referenceFile));
    if(!refFile.good()) {
        std::cerr << "ERROR: cannot open file of references\n";
        return 1;
    }
    std::vector<std::string> references;
    for( std::string line; std::getline( refFile, line ); )
    {
        references.push_back(line);
    }

    Graph g;
    std::vector<std::string>::iterator iter;
    for(iter = references.begin(); iter != references.end(); iter++) {

        // Translate from reference name to rID.
        int rID = 0;
        if (!getIdByName(rID, seqan::contigNamesCache(seqan::context(input_bam)), *iter))
        {
            std::cerr << "ERROR: Reference sequence named " << *iter << " not known.\n";
            return 1;
        }

        // Translate BEGIN and END arguments to number, 1-based to 0-based.
        int beginPos = 0, endPos = options.endLength - 1;
        if(parseRegion(rID, beginPos, endPos, bamContext, baiIndex, input_bam, g ) != 0)
        {
            return 1;
        }
        beginPos = seqan::contigLengths(bamContext)[rID] - options.endLength;
        endPos   = seqan::contigLengths(bamContext)[rID];
        if(parseRegion(rID, beginPos, endPos, bamContext, baiIndex, input_bam, g) != 0)
        {
            return 1;
        }
    
    }
    g.removeLinks(options.lower, options.upper);
    std::cout << g;
    return 0;
}
