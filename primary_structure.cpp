#include <cassert>
#include <iostream>
#include "primary_structure.h"
#include "mesh.h"


int PrimaryStructure::newVertIndex(UvPos uv){
    uvMap.push_back(uv);
    return nv++;
}


bool operator == (const Half& a, const Half &b){ return a.next == b.next; }
bool operator != (const Half& a, const Half &b){ return !(a.next == b.next); }

int PrimaryStructure::maxVi() const{
    int res = -1;
    for (const AminoAcid &a:sequence) {
        res = std::max( res, a.vi );
    }
    return res;
}

void PrimaryStructure::clear(){
    nv = 0;
    sequence.clear();
    sequenceStart = 0;
}

PrimaryStructure reversePrefolding(const Mesh& m){
    PrimaryStructure res;

    int startHi = m.findStartingCutHalf(); // arbitrary!
    int hi = startHi;

    do {
        int valency;
        hi = m.followCut( hi , valency );
        res.sequence.push_back( AminoAcid(valency, hi ) );
    }while ( hi != startHi);

    res.initAllNextAndPrev();
    res.initVi(m);
    res.initAllPosAndDir();

    return res;
}

void PrimaryStructure::initAllNextAndPrev(){
    int n = sequence.size();
    for (int i=0;i<sequence.size();i++){
        sequence[i].next = (i+1)%n;
        sequence[i].prev = (i+n-1)%n;
    }
}

void PrimaryStructure::initVi(const Mesh& m){
    nv = 0;
    uvMap.clear();
    origToNewVi.clear();
    origToNewVi.resize(m.verts.size(),-1);
    for (AminoAcid &a:sequence) {

        int viOnMesh = m.halfs[a.hiOnMesh].vi;
        int vi = origToNewVi[viOnMesh];
        if (vi == -1) {
            vi = newVertIndex(a.pos());
            origToNewVi[viOnMesh] = vi;
        }
        a.vi = vi;
    }
}


void PrimaryStructure::initVi(){
    nv = 0;

    for (AminoAcid &a:sequence) {
        nv = std::max(nv,a.vi);
    }
    nv++;
}

void PrimaryStructure::splitAllCuts(){
    nv = 0;
    uvMap.clear();
    for (AminoAcid &a:sequence) {
        a.vi = newVertIndex(a.pos());
    }
}

void PrimaryStructure::initAllPosAndDir(){
    PrimaryStructure res;

    Tortoise startTortoise { UvPos(0,0), E }; // it's arbitrary!

    int i = sequenceStart;
    Tortoise tortoise = startTortoise;
    do {
        AminoAcid& a = sequence[i];
        a.tortoise = tortoise;
        tortoise.turn( a.val );
        tortoise.stepForward();
        i = a.next;
    }while ( i != sequenceStart);

    assert( (tortoise == startTortoise) && "Sequence not a closed 2D cycle!");

}

int PrimaryStructure::sequenceLength() const{
    int res = 0;
    int i = sequenceStart;
    do {
        i = sequence[i].next;
        res++;
        if (res>sequence.size()) return -1;
    } while (i!=sequenceStart);
    return res;
}

char valencyToChar(int val);

std::string PrimaryStructure::toString(int maxlen) const{
    std::string s;
    int i = sequenceStart;

    //s += std::to_string(sequenceStart)  + ",";
    //s += std::to_string(nv)  + " ";
    do {
        const AminoAcid& a = sequence[i];
        s += std::string() +  valencyToChar(a.val); //  + "" + std::to_string(a.vi) + "-";
        i = a.next;
        if (--maxlen==0) {
            s += std::string("...(") + std::to_string(sequenceLength()) + " aminos)";
            break;
        }
    } while (i!=sequenceStart);

    /*s+= "(";
    for (int i:origToNewVi) s += std::to_string(i) + " ";
    s+= ")";*/
    return s;
}

void PrimaryStructure::renameVerticesAs(const PrimaryStructure& other){
    int s = sequence.size();
    assert(s==other.sequence.size() );

    std::vector<int>  newNames(s,-1);

    for (int i=0; i<s; i++) {

        int tvi = sequence[i].vi;
        int tvj = other.sequence[i].vi;

        sequence[i].vi =  other.sequence[i].vi;

        if (newNames[tvi] == -1) newNames[tvi] = tvj;
        else assert(newNames[tvi] == tvj);
    }
    for (int &i:origToNewVi){
        if (i>=0) {
            assert(newNames[i]>=0);
            i = newNames[i];
        }
    }

}

std::vector<int> PrimaryStructure::vertexMapTo(const PrimaryStructure& other) const{
    int s = sequence.size();
    assert(s==other.sequence.size() );
    assert(sequenceStart == other.sequenceStart );
    std::vector<int> res( s , -1 );
    int maxi=0, maxj=0;
    for (int i=0; i<s; i++) {
        int vi = sequence[i].vi;
        int vj = other.sequence[i].vi;
        maxi = std::max(vi,maxi);
        maxj = std::max(vj,maxj);
        if (res[vi] == -1) res[vi] = vj;
        else assert(res[vi] == vj);
    }
    assert(maxi==maxj);
    res.resize(maxi+1);
    return res;
}

