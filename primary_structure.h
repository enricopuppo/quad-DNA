#pragma once
#include<vector>
#include<string>
#include "mesh.h"
#include "uv.h"



struct Tortoise{
    UvPos p;
    UvDir d;
    void stepForward(int nSteps){
        for (int i=0; i<nSteps; i++) stepForward();
    }
    void stepForward(){
        p = p+d;
    }
    void stepBackward(){
        p = p-d;
    }
    void turn(int valnecy) { d = d+valnecy;}

    bool operator ==(const Tortoise &other) const { return p==other.p && d==other.d; }
    bool operator <(const Tortoise &other) const {
        return p.u<other.p.u ||
               (p.u==other.p.u && p.v<other.p.v ) ||
               (p.u==other.p.u && p.v==other.p.v && d<other.d );

    }
};


struct AminoAcid{
    bool isConcave() const { return val > 2; }
    bool isStraight() const { return val == 2; }
    bool isConvex() const { return val == 1; }

    int next,prev;  // links for a double-linked list (indices)

    Tortoise tortoise; // pos+dir of this boundary edge. Sits at the start.
    UvDir dir() const {return tortoise.d;}
    UvPos pos() const {return tortoise.p;}

    int val; // (valency) always >0. T A C G or beyond
    int vi; // vertex index. Arbitrary, but same index = same vertex.  (no gaps in numeration)

    int hiOnMesh; // half-edge index on mesh (only used in Primary Structure encoded from a mesh)
    AminoAcid(int val, int id):val(val),vi(id),hiOnMesh(id){}

};

struct Quad;
class RNA;

class PrimaryStructure
{
public:

    bool exportAsObj(const std::string& objFilename) const;
    bool exportUvMapAsObj(const std::vector<Quad>& quads, const std::string& objFilename) const;

    std::string toString(int maxlen = 0) const;

    std::vector<Quad> prefolding(const Mesh* m = NULL);

    // data linking to orig mesh
    std::vector<int> origToNewVi;
    std::vector<UvPos> uvMap; // per created vert

    void splitAllCuts(); // call me to build a good uvMap

    std::vector<int> vertexMapTo(const PrimaryStructure& other) const;
    void renameVerticesAs(const PrimaryStructure& other);


private:

    int sequenceStart = 0;

    std::vector<AminoAcid> sequence;
    int nv = 0;

    void initAllNextAndPrev();
    void initVi(const Mesh& m);
    void initVi();
    void initAllPosAndDir();

    int newVertIndex(UvPos uv);

    // local moves for secondary structure construction
    Quad bonkRetract(AminoAcid &a,AminoAcid &b,AminoAcid &c,AminoAcid &d);
    Quad cornerRetract(AminoAcid &a,AminoAcid &b,AminoAcid &c);
    Quad finalQuad(AminoAcid &a,AminoAcid &b,AminoAcid &c,AminoAcid &d) const;
    void decreaseValency(AminoAcid &a);
    void assignValency(AminoAcid &a, int newVal);
    void removePairFromSequence(AminoAcid &a, AminoAcid &b);
    void setInitialRedFlags();

    bool isRedFlagged(const AminoAcid &a) const;
    void addRedFlags(const AminoAcid &a);
    void removeRedFlags(const AminoAcid &a);

    int maxVi() const;

    void clear();

    int sequenceLength() const;

friend class Mesh;

friend PrimaryStructure translation(const RNA &);
friend RNA reverseTranslation(const PrimaryStructure &);

friend PrimaryStructure reversePrefolding(const Mesh&);
};

PrimaryStructure reversePrefolding(const Mesh& m);

