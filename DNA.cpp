#include <sstream>
#include <iostream>
#include <assert.h>
#include "DNA.h"
#include "RNA.h"

char valencyToChar(int val);

std::string DNA::toString() const{
    std::ostringstream s;
    for (const BasePair& b:helix){
        s   << valencyToChar(b.base[0]) << valencyToChar(b.base[1])  << b.len;
    }
    return s.str();
}



static std::string separatorString(const BasePair&p){
    static const char* openPars[] = {
        "",
        "(",
        "((",
        "(((",
        "((((",
        "(((((",
        "((((((",
        "(((((((",
        "((((((((",
    };

    char sep = '|';
    int k = p.openRounds;
    if (p.closedRound) {
        if (p.portOut) sep = '}'; else sep = ')';
    } else {
        if (k>0) {
            sep = '(';
            k--;
        }
    }
    return std::string() + sep + openPars[k];
}

static std::string portsInString(const BasePair&p){
    std::string res;
    for (int portIn:p.portsIn) res += 'a' + portIn;
    return res;

}

template<class T>
std::vector<T> operator + (const std::vector<T>& a,const std::vector<T>& b){
    std::vector<T> c = a;
    for (const T& bi:b) c.push_back(bi);
    return c;
}

BasePair operator + (const BasePair& a,const BasePair& b){
    BasePair c;
    c.base[0] = a.base[0] + b.base[0] - 2;
    c.base[1] = a.base[1] + b.base[1] - 2;
    c.len = a.len + b.len;
    c.openRounds = a.openRounds + b.openRounds;
    c.closedRound = a.closedRound || b.closedRound;
    c.portOut = b.portOut || a.portOut;
    c.portsIn = b.portsIn + a.portsIn ;
    assert(a.len==0 || b.len==0 );
    assert(!a.closedRound || !b.closedRound);
    return c;
}

void operator += (BasePair& a,const BasePair& b){
    a = a+b;
}


bool BasePair::isCanon() const{
    if (openRounds>1) return false;
    if (portsIn.size()>1) return false;
    if (base[0]>4) return false;
    if (base[1]>4) return false;  // no baes beyond G ("go back")
    return true;
}

static void pop_front(std::vector<int> &v){
    for (int i=0; i<v.size()-1; i++) v[i] = v[i+1];
    v.pop_back();
}

BasePair BasePair::splitCanonPad(){
    BasePair res;
    res.openRounds = 0;
    res.closedRound = false;
    res.len = 0;
    res.portOut = false;

    if (openRounds>1) {
        res.openRounds = 1;
        openRounds--;
    }

    if (portsIn.size()>1) {
        res.portsIn.push_back( portsIn[0] );
        pop_front(portsIn);
    }
    for (int i=0; i<2; i++) if (base[i]>4) {
        int excessVal = std::min(base[i]-4,2);
        res.base[i] = 2+excessVal;
        base[i]-=excessVal;
    } else res.base[i] = 2;
    assert(res.isCanon());

    return res;
}


DNA DNA::canonicalVersion() const{
    DNA res;

    for (BasePair b:helix){
        while (!b.isCanon()) {
            res.helix.push_back(b.splitCanonPad());
        }
        res.helix.push_back(b);
    }
    return res;
}
void BasePair::setZero(){
    len = 0;
    closedRound = false;
    base[0] = base[1] = 2;
    openRounds = 0;
    portOut = 0;
    portsIn.clear();
}


DNA DNA::nonCanonicalVersion() const{

    DNA res;

    BasePair cumulated;
    cumulated.setZero();

    for (BasePair b:helix){
        if (b.len==0) {
            cumulated += b;
        } else {
            res.helix.push_back(b+cumulated);
            cumulated.setZero();
        }
    }

    return res;
}

bool DNA::fromString(const std::string&){
    return false;
}

std::string DNA::toStringFat() const{


    std::ostringstream s;
    const char* indent = "\n     ";
    //s << indent;

    for (const BasePair& b:helix){
        s   << valencyToChar(b.base[0]) << " ";
    }
    s << indent;
    for (const BasePair& b:helix){
        s  << separatorString(b) << portsInString(b) << std::to_string(b.len) ;
    }
    s << indent;
    for (const BasePair& b:helix){
        s   << valencyToChar(b.base[1]) << " ";
    }    return s.str();
}


std::string DNA::toStringArt(const char* prefix) const{
    std::vector<int> basePos(helix.size(),0);
    int asciiLen = 0;
    int k=0;
    std::ostringstream midS;
    midS << prefix;
    asciiLen += std::string(prefix).length();
    for (const BasePair& b:helix) {
        std::string sep = separatorString(b);
        std::string len = portsInString(b) + std::to_string(b.len);
        basePos[k++] = asciiLen + sep.length()/2;
        asciiLen += sep.length()+len.length();
        midS << sep << len;
    }
    if (helix[0].portOut) midS << "}";

    midS << '\n';
    asciiLen++;

    std::string above(asciiLen,' ');
    std::string below(asciiLen,' ');
    above[asciiLen-1] = '\n';
    below[asciiLen-1] = '\n';

    k = 0;
    for (const BasePair& b:helix) {
        int p = basePos[k++];
        above[p] = valencyToChar(b.base[0]);
        below[p] = valencyToChar(b.base[1]);
    }

    std::ostringstream s;
    s << above << midS.str() << below;
    return s.str();
}

/* Auxiliary data for reverseTrascritions */

struct VisitedNode{
    int spanId, nodeId, helixPos;
    VisitedNode(Chunk c):spanId(c.spanId),nodeId(c.nodeId),helixPos(-2) {};
};

static std::vector<VisitedNode> visited;

static bool alreadyVisited(int nodeId){
    if (nodeId == 0) return true;
    for (const VisitedNode &v:visited) if (v.nodeId == nodeId) return true;
    return false;
}

static std::vector<int> flipAndClear(std::vector<int> &v){
    std::vector<int> res;
    res.reserve(v.size());
    for (int i=v.size()-1; i>=0; i--) res.push_back(v[i]);
    v.clear();
    return res;
}
/* end */

DNA reverseTranscription(const RNA & rna){
    DNA dna;

    int cn = rna.chunks.size();

    // for each span id: is this spanId connected to a port? if not: -1. If so: port number
    std::vector<int> openPorts(rna.maxSpanId(),-1);

    int genus = rna.evaluateGenus();
    if (genus != 0) std::cout << "DNA for Genus "<< genus << " >0 ! Let's go!!!" << std::endl;

    visited.clear();

    enum{DOWN, UP} direction = DOWN;

    int baseOfLastLeaf = -1;
    bool lastLeafWasTakeOff = false;
    std::vector<int> lastPortIns;
    int nPorts = 0;

    int nIte = 0;
    int safety = 500;

    for (int ci=0; ; ci=(ci+1)%cn) {
        if (safety--<0) break;
        const Chunk& c = rna.chunks[ci];
        int destNodeId = rna.chunks[ (ci+1) % cn].nodeId;


        if (!visited.empty() && visited.back().spanId == c.spanId) {
            // go up on the tree!
            //std::cout << "UP [" << ci << "] spanId=" <<  c.spanId <<std::endl;

            assert(visited.back().nodeId == destNodeId );

            if (nIte>=rna.chunks.size()) break;

            if (direction == DOWN) {
                // DOWN => UP: a leaf
                baseOfLastLeaf = c.nodeVal;
            } else {
                // UP => UP
                BasePair& b = dna.helix[visited.back().helixPos];
                assert( b.base[1] == -1);
                b.portsIn = flipAndClear(lastPortIns); //lastPortIns.clear();
                b.base[1] = c.nodeVal;
            }

            visited.pop_back();

            direction = UP;
        } else {
            // go down the tree!
            //std::cout << "DOWN -- rna[" << ci << "].spanId=" <<  c.spanId <<std::endl;


            // maybe we encountered an open port?
            if (openPorts[c.spanId]>=0 ) {
                ci = rna.otherIndexOfSameSpan(ci);
                //std::cout << "   !!!SKIP  ==> [" << ci << "] spanId=" << rna.chunks[ci].spanId <<std::endl ;
                lastPortIns.push_back( openPorts[c.spanId] );
                //if (lastPortIns.size()==2) break;// MEGA HACK JUST FOR NOW:
                openPorts[c.spanId] = -1;
                continue;
            }

            if (nIte>=rna.chunks.size()) break;

            BasePair b;
            b.portsIn = flipAndClear(lastPortIns); //lastPortIns.clear();
            b.base[0] = c.nodeVal;
            if (direction == DOWN) {
                // DOWN => DOWN
                b.base[1] = -1; // to be filled on the way back!
                if (!visited.empty()) visited.back().helixPos = dna.helix.size();
            } else {
                // UP => DOWN (an "armpit")
                if (!visited.empty()) dna.helix[visited.back().helixPos].openRounds++;
                b.base[1] = baseOfLastLeaf;
                b.closedRound = true;
                b.portOut = lastLeafWasTakeOff;
                lastLeafWasTakeOff = false;
            }
            if (alreadyVisited( destNodeId ) ) {
                lastLeafWasTakeOff = true;
                openPorts[ c.spanId ] = nPorts++;
                ci = rna.otherIndexOfSameSpan( ci );
                //std::cout << ">>>>>>Cartello (node " << destNodeId << " already vis) ";
                std::cout << " rna["<<ci<<"].spanId = " << rna.chunks[ci].spanId << std::endl;
                ci--;
            }
            b.len = rna.edgeLenghtOf(c);
            visited.push_back(c);
            dna.helix.push_back( b );
            direction = DOWN;
        }
        nIte++;
    }

    // is there are leftover in lastPortsIns, attach them to start
    BasePair& b = dna.helix[0];
    b.portsIn = flipAndClear(lastPortIns); //lastPortIns.clear();


    std::cout << nPorts << " teleports (genus " << nPorts/2 << "?)" << std::endl;

    assert( dna.helix[0].base[1] == -1);
    dna.helix[0].base[1] = baseOfLastLeaf;
    dna.helix[0].portOut = lastLeafWasTakeOff;

    std::cout<<dna.toStringArt("XXX =");

    return dna.canonicalVersion();
    //return dna;

}

/*
 * Simpler version, for Genus 0 only
 *
DNA reverseTranscription(const RNA & rna){
    DNA dna;

    if (rna.evaluateGenus() != 0) {
        std::cerr << "DNA not done. For now, I only can work with Genus 0 :-(" << std::endl;
        return dna;
    }

    visited.clear();
    enum{DOWN, UP} direction = DOWN;
    int baseOfLastLeaf = -1;

    for (int ci=0; ci<rna.chunks.size(); ci++) {
        const Chunk& c = rna.chunks[ci % rna.chunks.size()];
        const Chunk& cNext = rna.chunks[ (ci+1) % rna.chunks.size()];

        if (!visited.empty() && visited.back().spanId == c.spanId) {
            // go up on the tree!
            assert(visited.back().nodeId == cNext.nodeId );

            if (direction == DOWN) {
                // DOWN => UP: a leaf
                baseOfLastLeaf = c.nodeVal;
            } else {
                // UP => UP
                //std::cout << visited.back().helixPos << std::endl;
                BasePair& b = dna.helix[visited.back().helixPos];
                assert( b.base[1] == -1);
                b.base[1] = c.nodeVal;
            }
            visited.pop_back();

            direction = UP;
        } else {
            // go down the tree!

            BasePair b;
            b.base[0] = c.nodeVal;
            if (direction == DOWN) {
                // DOWN => DOWN
                b.base[1] = -1; // to be filled on the way back
                if (!visited.empty()) visited.back().helixPos = dna.helix.size();
            } else {
                // UP => DOWN (an "armpit")
                if (!visited.empty()) dna.helix[visited.back().helixPos].openRound++;
                b.base[1] = baseOfLastLeaf;
                b.closedRound = true;
            }
            b.len = rna.edgeLenghtOf(c);
            visited.push_back(c);
            dna.helix.push_back( b );

            direction = DOWN;
        }
    }

    assert( dna.helix[0].base[1] == -1);
    dna.helix[0].base[1] = baseOfLastLeaf;


    return dna;

}*/

DNA::Site DNA::step(DNA::Site s) const{
    if (s.side==0) return stepForward(s);
    else return stepBackward(s);
}

std::string DNA::Site::toString() const{
    return std::string(" (") + std::to_string(index+1) + ")" + std::string((side==0)?"=>":"<=");
}

std::pair<int, int> DNA::ithInPort(int i) const{
    for (int j=0; j<helix.size(); j++) {
        for (int k=0; k<helix[j].portsIn.size(); k++)
            if (helix[j].portsIn[k]==i) return std::make_pair(j,k);
    }
    assert(false && "i-th in port not found");
    return std::make_pair(-1,-1);
}

int DNA::ithOutPort(int i) const{
    int count = 0;
    //std::cout << "seeking " << (i+1) << "th '}'..." << std::endl;
    for (int j=1; j<helix.size(); j++) {
        if (helix[j].portOut) {
            if (count==i) {
                //std::cout << "found [" << j << "]" << std::endl;
                return j;
            }
            count++;
        }
    }
    return 0; // it's the last }, at index 0
    //assert(false && "i-th out-port not found");
    //return -1;
}

int DNA::rankOfOutPort(int helixPos) const{
    assert(helix[helixPos].portOut);
    int res = 0;
    for (int j=1; j<helixPos; j++) {
        if (helix[j].portOut) res++;
    }
    return res;
}


DNA::Site DNA::exitOfInPort(Site s, int which) const{
    int portIn= helix[s.index].portsIn[which];
    s.index = ithOutPort( portIn );
    return s;
}

DNA::Site DNA::exitOfLastInPort(Site s) const{
    int portIn = helix[s.index].portsIn.back();
    s.index = ithOutPort( portIn );
    return s;
}


DNA::Site DNA::exitOfOutPort(Site s) const{
    std::pair<int,int> p = ithInPort( rankOfOutPort( s.index ) );
    int newPos = p.first;
    int k = p.second;
    s.index = newPos;
    if (k>0) {
        s.indexOfVi = s.index;
        s = exitOfInPort( s, k-1 ); // two (o more) consecutive in ports!
    } else {
        s.indexOfVi = s.index;
        if (helix[s.index].closedRound) {
            // this out port is in an "armpit"
            s.side = 0;
            findMatchingOpenParOf( s.indexOfVi );
        }
    }

    return s;
}


DNA::Site DNA::stepForward(DNA::Site s) const{
    s.index++;
    BasePair b = helix[s.index % helix.size()];

    s.indexOfVi = s.index;

    if (b.closedRound || s.index == helix.size() ) {
        s.side = 1; // turn back
        if (b.portOut) {
            s = exitOfOutPort(s);
        }
    }

    return s;
}

bool DNA::findMatchingOpenParOf(int &helixPos) const{
    // run backward to matching open par; skip any "(...)"
    assert( helix[helixPos].closedRound );
    int inside = 1;
    do  {
        helixPos--;
        if (helix[helixPos].closedRound) inside++;
        inside -= helix[helixPos].openRounds;
        if (helixPos==0) inside--; // implict parentesis at the end
    } while (inside>0);
    return (inside<0);
}

DNA::Site DNA::stepBackward(DNA::Site s) const{
    // march backward by 1, skipping any "()" pair

    s.index--;
    if (s.index<0) s.index = helix.size()-1;

    int oldPos = s.index;
    bool openPar;

    if (helix[s.index].closedRound) {
        openPar = findMatchingOpenParOf( s.index );
    } else {
        openPar = helix[s.index].openRounds>0 || (s.index==0);
    }

    s.indexOfVi = s.index;

    if (openPar) {
        // bounce back
        s.index = oldPos;
        s.side = 0;
        // run Forward To matching Closed Parentesis; skip any "(...)"
        int inside=1;
        do  {
            s.index++;
            if (helix[s.index].closedRound) inside--;
            inside += helix[s.index].openRounds;
            s.index = s.index % helix.size();
            if (s.index==0) break;
        } while (inside>0 );
        if (s.index==0) s.endReached = true;
    }

    if (!helix[s.index].portsIn.empty()) {
        s = exitOfLastInPort( s );
        s.side = 1;
    }
    return s;
}

/*
 * Simpler version for genus 0 only
DNA::Site DNA::stepForward(DNA::Site s) const{
    s.index++;
    BasePair b = helix[s.index % helix.size()];
    if (b.closedRound || s.index == helix.size() ) {
        // turn back
        s.side = 1;
    }
    s.indexOfVi = s.index;
    return s;
}

DNA::Site DNA::stepBackward(DNA::Site s) const{
    // march backward by 1, skipping any "()" pair

    s.index --;

    int oldPos = s.index;
    bool openPar;

    if (helix[s.index].closedRound) {
        // run backward to matching open par
        int inside = 1;
        do  {
            s.index--;
            if (helix[s.index].closedRound) inside++;
            inside -= helix[s.index].openRound;
            if (s.index==0) inside--;
        } while (inside>0);
        openPar = (inside<0);
    } else {
        openPar = helix[s.index].openRound>0 || (s.index==0);
    }

    s.indexOfVi = s.index;

    if (openPar) {
        s.index = oldPos;
        s.side = 0;
        // run Forward To matching Closed Parentesis;
        int inside=1;
        do  {
            s.index++;
            if (helix[s.index].closedRound) inside--;
            inside += helix[s.index].openRound;
            s.index = s.index % helix.size();
            if (s.index==0) break;
        } while (inside>0 );
        if (s.index==0) s.endReached = true;
    }
    return s;
}
*/

/* Auxiliary data for trascritions */
static std::vector<int> viRemap;


Chunk DNA::chunkAt(Site s) const{
    Chunk c;
    int helixSize = helix.size();
    c.nodeVal = helix[s.index%helixSize].base[ s.side ];
    c.nodeId = viRemap[ s.indexOfVi  ];
    //assert(c.nodeId!=-1 && "Why the index of a } is being used as vertex index??");
    c.spanId = (s.side==0)? s.index: ((s.index+helixSize-1)%helixSize);
    return c;
}

RNA transcription(const DNA &d){
    const DNA dna = d.nonCanonicalVersion();
    std::cout<<dna.toStringArt("XXX =");

    RNA rna;

    viRemap.resize(dna.helix.size()+1,0);
    for (int i=0,c=0; i<viRemap.size(); i++) viRemap[i] = (i!=0 && dna.helix[i%dna.helix.size()].portOut) ? (-1) : (c++);
    for (int H:viRemap) std::cout << " " << H; std::cout << std::endl;

    DNA::Site sStart({0,0,0,false});
    DNA::Site s = sStart;

    //int safety = 500;
    do {
        Chunk c = dna.chunkAt(s);
        rna.chunks.push_back( c );
        if (rna.chunks.size() == dna.helix.size()*2 ) break;
        s = dna.step(s);
        //if (safety--<0) break;
    } while (1);

    for (BasePair b : dna.helix) rna.spanLenghts.push_back( b.len );

    return rna;
}



