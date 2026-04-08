#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <glm/vec3.hpp>
#include "union_find.h"

// qi: quad index (a face of the mesh)
// vi: vertex index
// ni: node index (a node of the cut-graph)
// hi: half index (an half-edge)
// ei: edge index

class PrimaryStructure;
class UvPos;

struct Half
{
    int next;
    int opp;
    int qi;
    int ei;
    int vi; // emanating vert

    int mark;
};

struct Edge
{
    int hi;
    bool isCut;
    int mark;
};

struct Quad
{
    int vi[4];
    glm::vec3 norm;
    int hi;
};

struct Vert
{
    int hi;
    int valency;
    int ni = -1; // node index, -1 if not a node
    bool isRegular() const { return valency == 4; }
    bool isNode() const { return ni >= 0; }

    glm::vec3 pos;
    int mark;
};

// cut graph
struct Node
{
    int vi;
    int valency; // of arcs
};

struct Arc
{
    int ni, nj;          // node indices
    std::vector<int> ei; // edge indices
    float priority;
    bool dissolved;
    bool operator<(const Arc &other) const { return priority > other.priority; }
};

struct BoundingCube
{
    glm::vec3 center;
    float radius;
};

class Mesh
{
    std::vector<Half> halfs;

    std::vector<Quad> quads;
    std::vector<Edge> edges;
    std::vector<Vert> verts;

    // cut-graph
    std::vector<Node> nodes;
    std::vector<Arc> arcs;
    UnionFind quadRegions;

    // logging stream, defaults to std::cout. File stream lifetime is managed outside Mesh.
    std::ostream *logFile = &std::cout;

    // construction
    void populateHalfs();
    void populateVerts(); // not necessary if imported from obj
    void populateEdges();

    // update
    void updateVertexValencies();
    void clearCuts();

    // marking elements
    int incrementalMark = 0;
    template <class Type>
    bool isMarked(const Type &t) const { return t.mark == incrementalMark; }
    template <class Type>
    void mark(Type &t) { t.mark = incrementalMark; }
    void unmarkAll() { incrementalMark++; }

    // access
    Vert vert(Half h) const { return verts[h.vi]; }
    Vert &vert(Half h) { return verts[h.vi]; }
    Vert vertEnd(Half h) const { return verts[halfs[h.opp].vi]; }
    Vert &vertEnd(Half h) { return verts[halfs[h.opp].vi]; }
    int vertEndI(Half h) const { return halfs[h.opp].vi; }
    Quad quad(Half h) const { return quads[h.qi]; }
    Quad &quad(Half h) { return quads[h.qi]; }
    Edge edge(Half h) const { return edges[h.ei]; }
    Edge &edge(Half h) { return edges[h.ei]; }
    int viOf(const Edge &e) const { return halfs[e.hi].vi; }
    int vjOf(const Edge &e) const { return opposite(halfs[e.hi]).vi; }
    int qiOf(const Edge &e) const { return halfs[e.hi].qi; }
    int qjOf(const Edge &e) const { return opposite(halfs[e.hi]).qi; }

    // navigation
    int opposite(int hi) const { return halfs[hi].opp; }
    Half opposite(const Half &h) const { return halfs[h.opp]; }
    Half &opposite(const Half &h) { return halfs[h.opp]; }
    int next(int hi) const { return halfs[hi].next; }
    Half next(const Half &h) const { return halfs[h.next]; }
    Half &next(const Half &h) { return halfs[h.next]; }
    int pivotCW(int hi) const { return next(opposite(hi)); }
    Half pivotCW(const Half &h) const { return next(opposite(h)); }
    Half &pivotCW(const Half &h) { return next(opposite(h)); }
    int ahead(int hi) const { return next(opposite(next(hi))); }
    Half &ahead(const Half &h) { return next(opposite(next(h))); }

    // methods used by cut-graph optimization
    void populateNodes();
    void populateArches();
    bool dissolveArches();
    void assignQuadsToRegions();
    float evaluatePriority(Arc &a) const;
    // auxiliary:
    void createArchFrom(Half h);
    bool canBeDissolved(Arc a);
    void dissolve(Arc &a);
    bool compressArcsVector();

    // methods used by flattening:
    int findStartingCutHalf() const;
    int followCut(int hi, int &valency) const;

    // diagnostics
    void print(const Arc &a) const;
    void print(const Node &n) const;
    void printNodes() const;

    // attibute vectors
    void scatterVertAttr(const std::vector<float> &perFace, std::vector<float> &perVert) const;
    void gatherQuadAttr(const std::vector<float> &perVert, std::vector<float> &perQuad) const;

    // functionalities
    friend class MCG;
    friend class GlViewer;
    friend class PrimaryStructure;
    friend PrimaryStructure reversePrefolding(const Mesh &m);

public:
    void setLogFile(std::ostream &out) { logFile = &out; }
    std::ostream &getLogFile() { return *logFile; }
    const std::ostream &getLogFile() const { return *logFile; }

    bool importObj(const std::string &filename);
    bool exportObj(const std::string &filename) const;
    bool exportObj(const std::string &filename, const Mesh &uvMap) const;
    bool exportCutsAsObj(const std::string &filename) const;

    void computeMotorcycleGraphCuts();
    void dissolveCutsToMakeOneDisk();

    void assignConnectivity(const std::vector<Quad> &q) { quads = q; }
    void assignEmbeddingAs(const Mesh &other);
    void assignEmbeddig(const std::vector<UvPos> &uvmap);
    void normalizeIn01();

    void renameVertices(const std::vector<int> newVi);

    // stats & diagnostics
    int nRegions() const { return quadRegions.nClusters(); }
    bool isClosed() const;
    int countIrregulars() const;
    BoundingCube computeBoundingCube() const;

    std::string statsString() const;
};
