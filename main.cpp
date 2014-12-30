#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include "version.h"

// graph typedefs - the cargo is information for the edges

// typedefs for the bam files
typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

struct Options {
    int upper;
    int lower;
    seqan::CharString bamFile;
    seqan::CharString baiIndexFile;
    seqan::CharString referenceFile;

    Options() : upper(-1), lower(3)
    {}
};

class Graph {
    public:
        void addLink(seqan::String<char> id1, seqan::String<char> id2);
        void removeLinks(int lower = 3, int upper = -1);
    private:
        std::map<std::pair<seqan::String<char>, seqan::String<char> >, int > adjacency_list;

        friend std::ostream& operator<<(std::ostream& os, Graph& g);

};
void Graph::addLink(seqan::String<char> id1, seqan::String<char> id2) {
    std::pair<seqan::String<char>, seqan::String<char> > key;
    if( id1 < id2) {
        key = std::make_pair(id1,id2);
    } else {
        key = std::make_pair(id2,id1);
    }
    if(adjacency_list.find(key) == adjacency_list.end()) {
        adjacency_list[key] = 1;
    } else {
        adjacency_list[key] += 1;
    }
}

void Graph::removeLinks(int lower, int upper) {
    std::map<std::pair<seqan::String<char>, seqan::String<char> >, int >::iterator iter;
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
    std::map<std::pair<seqan::String<char>, seqan::String<char> >, int >::iterator iter;
    for(iter = g.adjacency_list.begin(); iter != g.adjacency_list.end(); ++iter) {
        os <<iter->first.first<<"\t"<<iter->first.second<<"\t"<<iter->second<<std::endl;
    }
    return os;
}

int parseRegion(int rID, int beginPos, int endPos, TBamIOContext& context, seqan::BamIndex<seqan::Bai>& baiIndex, seqan::BamHeader& header, seqan::Stream<seqan::Bgzf>& inStream, Graph& g ) {
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inStream, hasAlignments, context, rID, beginPos, endPos, baiIndex))
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
        if (readRecord(record, context, inStream, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM file.\n";
            return 1;
        }

        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || record.rID != rID || record.beginPos >= endPos)
            break;
        // If we are left of the selected position then we skip this record.
        if (record.beginPos < beginPos)
            continue;

        // make sure that the mate is not on the contig
        if(record.rID != record.rNextId && record.rNextId != -1) {
            // Otherwise, we add it to the graph.
            g.addLink(header.sequenceInfos[record.rID].i1, header.sequenceInfos[record.rNextId].i1);
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
    seqan::addOption(parser, seqan::ArgParseOption( "b", "bam", "Bam file", seqan::ArgParseArgument::INPUTFILE, "FILE"));
    seqan::setRequired(parser, "bam");
    seqan::addOption(parser, seqan::ArgParseOption( "B", "bai", "Index file", seqan::ArgParseArgument::INPUTFILE, "FILE"));
    seqan::setRequired(parser, "bai");
    seqan::addOption(parser, seqan::ArgParseOption( "r", "ref-seqs", "Contigs to look for", seqan::ArgParseArgument::INPUTFILE, "FILE"));
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

    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, seqan::toCString(options.bamFile), "r"))
    {
        std::cerr << "ERROR: Could not open " << options.bamFile << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (read(baiIndex, seqan::toCString(options.baiIndexFile)) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << options.baiIndexFile << "\n";
        return 1;
    }

    // Setup name store, cache, and BAM I/O context.
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, inStream, seqan::Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file " << options.bamFile << "\n";
        return 1;
    }

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
        if (!getIdByName(nameStore, *iter, rID, nameStoreCache))
        {
            std::cerr << "ERROR: Reference sequence named " << *iter << " not known.\n";
            return 1;
        }

        // Translate BEGIN and END arguments to number, 1-based to 0-based.
        int beginPos = 0, endPos = 499;
        if(parseRegion(rID, beginPos, endPos, context, baiIndex, header, inStream, g ) != 0)
        {
            return 1;
        }
        beginPos = header.sequenceInfos[rID].i2 - 500;
        endPos   = header.sequenceInfos[rID].i2;
        if(parseRegion(rID, beginPos, endPos, context, baiIndex, header, inStream, g) != 0)
        {
            return 1;
        }
    
    }
    g.removeLinks(options.lower, options.upper);
    std::cout << g;
    return 0;
}
