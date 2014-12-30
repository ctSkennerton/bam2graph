#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

// graph typedefs - the cargo is information for the edges

// typedefs for the bam files
typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

class Graph {
    public:
        void addLink(seqan::String<char> id1, seqan::String<char> id2);
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


int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam IN.bam.bai Ref_file\n";
        return 1;
    }

    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, argv[1], "r"))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (read(baiIndex, argv[2]) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << argv[2] << "\n";
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
        std::cerr << "ERROR: Could not read header from BAM file " << argv[1] << "\n";
        return 1;
    }

    // read the references file
    std::ifstream refFile(argv[3]);
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
            std::cerr << "ERROR: Reference sequence named " << argv[3] << " not known.\n";
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
    std::cout << g;
    return 0;
}
