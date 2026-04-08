#include <sstream>
#include <iostream>
#include <map>
#include "RNA.h"
#include "primary_structure.h"


char valencyToChar(int val){
    switch (val) {
    case 1: return 'T'; // Turn (convex)
    case 2: return 'A'; // Ahead
    case 3: return 'C'; // Concave turn
    case 4: return 'G'; // Go back
    case 5: return 'U'; // extreme convex valency...
    case 6: return 'V';
    case 7: return 'W';
    case 8: return 'X';
    case 9: return 'Y';
    case 10: return 'Z';
    default: {
        // just for marked RNA/DNA (for debug)
        if (val<1) return '?';
        const int nLetters = 26;
        int letter = (val-11)%(nLetters*2);
        if (letter<nLetters) return 'A'+letter;
        else if (letter<nLetters*2) return 'a'+(letter-nLetters);
        else return '!'; //132 + (letter-nLetters*2);
    }
    }
}

int charToValency(char c){
    switch (c) {
    case 'A': return 2; // Ahead
    case 'C': return 3; // Concave turn
    case 'G': return 4; // Go back
    case 'T': return 1; // Turn (convex)
    case 'U': return 5; // extreme convex valencies
    case 'V': return 6;
    case 'W': return 7;
    case 'X': return 8;
    case 'Y': return 9;
    case 'Z': return 10;
    default: return 0;
    }
}

std::string RNA::toString() const{

    std::ostringstream s;
    const char* endingStr = "";
    for (const Chunk& c:chunks){
        s  << endingStr << c.nodeId << valencyToChar(c.nodeVal) << c.spanId;
        endingStr = "|";
    }

    s << "[";
    endingStr = "";
    for (int len:spanLenghts) {
        s  << endingStr << len;
        endingStr = ",";
    }
    s << "]";
    return s.str();
}

bool RNA::fromString(const std::string& st){
    clear();
    std::istringstream s(st);
    while (!s.eof()) {
        char sep, cVal;
        Chunk c;
        s >> c.nodeId >> cVal >> c.spanId >> sep ;
        c.nodeVal = charToValency(cVal);
        chunks.push_back( c );
        if (sep == '[') break;
        if (sep != '|') goto wrong;
    }
    while (!s.eof()) {
        char sep;
        int i;
        s >> i >> sep;
        spanLenghts.push_back(i);
        if (sep == ']') break;
        if (sep != ',') goto wrong;
    }
    return true;
wrong:
    std::cerr << "Error loading RNA from \"" << st << "\"\n";
    return false;
}

int RNA::edgeLenghtOf(const Chunk &c) const{
    return vertLenghtOf(c)+1;
}

int RNA::vertLenghtOf(const Chunk &c) const{
    return spanLenghts[c.spanId]-1;
}

/*
 * translation:
 * PrimaryStructure <=> RNA
 */

PrimaryStructure translation(const RNA &rna) {
    PrimaryStructure p;
    p.clear();
    int ns = rna.spanLenghts.size();
    int firstAvailableNodeId = 0;
    for (Chunk c:rna.chunks){
        firstAvailableNodeId = std::max(c.nodeId+1, firstAvailableNodeId);
    }

    std::vector<int> startAt(ns,-1);
    for (Chunk c:rna.chunks){
        int si = c.spanId;
        int len = rna.vertLenghtOf(c);

        if (len == -1){
            // rotate in place
            assert(!p.sequence.empty() && "RNA shouldn't start with a zero-len span!");
            //p.sequence.back().val +=
        } else {

            p.sequence.push_back( AminoAcid( c.nodeVal, c.nodeId) );

            if (startAt[si]==-1) {
                // first encounter of this span
                startAt[si]=firstAvailableNodeId;
                for (int i=0; i<len; i++) {
                    p.sequence.push_back( AminoAcid(2,firstAvailableNodeId+i) );
                }
                firstAvailableNodeId += len;
            } else {
                // second (and last) encounter
                for (int i=0; i<len; i++) {
                    p.sequence.push_back( AminoAcid(2,startAt[si]+len-1-i) );
                }
            }
        }
    }

    p.uvMap.clear();
    p.origToNewVi.clear();

    p.initAllNextAndPrev();
    p.initVi();
    p.initAllPosAndDir();

    return p;
}

void RNA::clear(){
    spanLenghts.clear();
    chunks.clear();
}

RNA reverseTranslation(const PrimaryStructure & p){
    assert(p.sequenceStart==0 && p.sequence[0].next == 1 && "please call me before contracting me");
    assert(!p.sequence.empty());

    std::map<int,int> lastAminoacidToSpanId;

    int nv = p.maxVi()+1; // number of DISTINCT vertices

    // identify nodes (any turn, plus vertices not repeated exactly twice)
    std::vector<bool> isNode( nv, false);

    {
        std::vector<int> occurrences( nv, 0);
        for (AminoAcid a:p.sequence) {
            if (a.val!=2) isNode[a.vi] = true;
            occurrences[a.vi]++;
        }
        for (int i=0; i<nv; i++) if (occurrences[i]!=2) isNode[i] = true;
    }
    assert(isNode[0] && "primary structure should start from a node!");

    RNA rna;

    std::vector<int> viOfOrig(nv);
    std::vector<int> origOfVi;

    for (int viOrig=0; viOrig<nv; viOrig++) {
        if (isNode[viOrig]) {
            int viNew = origOfVi.size();
            origOfVi.push_back(viOrig);
            viOfOrig[viOrig] = viNew;
        }
    }

    std::vector<int> chunkStart;

    for (int i=0; i<p.sequence.size(); i++) {
        int origVi =  p.sequence[i].vi;
        if (isNode[origVi]) {
            chunkStart.push_back(i);
        }
    }

    int nChunks = chunkStart.size();

    rna.chunks.resize( nChunks );

    std::map<int,int> map; // from id of LAST vertex, to its span id


    for (int ci=0; ci<nChunks; ci++) {
        int si = chunkStart[ci]; // index in sequence
        const AminoAcid& a = p.sequence[ si ];
        rna.chunks[ci].nodeId = viOfOrig[a.vi];
        rna.chunks[ci].nodeVal = a.val;
        int siNext = (ci==rna.chunks.size()-1)? p.sequence.size() : chunkStart[ci+1];
        int len = siNext - si;
        int spanFirst = p.sequence[si+1].vi;
        int spanLast = p.sequence[siNext-1].vi;
        //assert(!isNode[spanFirst]);
        //assert(!isNode[spanLast]);

        auto found = map.find(spanFirst);
        if (found == map.end()) {
            // not found
            int newSpan = rna.spanLenghts.size();
            map[spanLast] = newSpan;
            rna.chunks[ci].spanId = newSpan;
            rna.spanLenghts.push_back(len);
        } else {
            rna.chunks[ci].spanId = found->second;
            assert(len == rna.spanLenghts[rna.chunks[ci].spanId]);
        }


    }


    return rna;
}


int RNA::evaluateGenus() const{
    int excessVal = 0;
    for (const Chunk& c:chunks) {
        excessVal += c.nodeVal - 2;
    }
    if( excessVal%4 != 0) {
        std::cerr << "Not a valid genus!" << std::endl;
        return -100;
    }
    return 1 + excessVal/4;
}

int RNA::otherIndexOfSameSpan(int i) const{
    int spanId = chunks[i].spanId;
    for (int j=0; j<chunks.size(); j++) {
        if (j!=i && chunks[j].spanId == spanId) return j;
    }
    assert(0 && "cannot find other instance with same Span Id (opposite)");
    return -1;
}

void RNA::markForDebug(){
    int k = 11;
    for (Chunk& c:chunks) {
        //std::cout << k << valencyToChar(k);
        c.nodeVal = k++;
    }

    for (int i=0; i<spanLenghts.size(); i++ ) {
        spanLenghts[i] = i;
    }
    std::cout << std::endl;
}



