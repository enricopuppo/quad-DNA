#include <algorithm>
#include <iostream>
#include "mesh.h"

void Mesh::dissolveCutsToMakeOneDisk(){
    std::cout << "\nSimplify graphs..." << std::endl;

    assignQuadsToRegions();
    populateNodes();
    populateArches();
    //printNodes();

    while (dissolveArches());

}

float Mesh::evaluatePriority(Arc& a) const{
    int vi = verts[nodes[a.ni].vi].valency;
    int vj = verts[nodes[a.nj].vi].valency;
    return ((vi==4) + (vj==4)) * 100.0 + // we prefer to dissolve arches stemming from regular vertices
           ((vi<4) + (vj<4)) * 10.0 +    // we prefer to keep arches stemming from high valency vertices
           a.ei.size() * 1.0 +           // we prefer to dissolve long arches
           a.nj*0.001 +                  // arbitary, to make priority not == and algo deterministic across platforms
           a.ni*0.011;
}

bool Mesh::canBeDissolved(Arc a){
    // (this method should be const but cannot be because of the union-find)

    // print(a);

    int arcValenceI = nodes[a.ni].valency;
    int arcValenceJ = nodes[a.nj].valency;
    bool regularI = verts[nodes[a.ni].vi].isRegular();
    bool regularJ = verts[nodes[a.nj].vi].isRegular();
    bool deadEndI = (arcValenceI == 1);
    bool deadEndJ = (arcValenceJ == 1);

    assert( arcValenceI>0 && arcValenceJ>0 );

    // rule A: at least one arc must be kept over any irregular vertex
    if (deadEndI && !regularI) {
        //std::cout<<" disallowed for rule A\n";
        return false;
    }
    if (deadEndJ && !regularJ) {
        //std::cout<<" disallowed for rule A\n";
        return false;
    }

    // rule B: cannot delete last arc
    if (deadEndI && deadEndJ) {
        //std::cout<<" disallowed for rule C\n";
        return false;
    }

    // rule C: no loops of quads, cut surface must be kept a disk!
    Edge e = edges[ a.ei[0] ]; // an arbitrary edge of the arc
    bool sameRegion = quadRegions.areSameCluster(qiOf(e),qjOf(e));
    if (!deadEndI && !deadEndJ && sameRegion) {
        //std::cout<<" disallowed for rule C\n";
        return false;
    }

    //std::cout<<" ARC IS KILLABLE!\n";
    return true;
}

void Mesh::dissolve(Arc &a){
    nodes[a.ni].valency--;
    nodes[a.nj].valency--;
    for (int i:a.ei) {
        assert(edges[i].isCut);
        edges[i].isCut = false;
    }
    Edge e = edges[a.ei[0]];
    quadRegions.fuse(qiOf(e),qjOf(e));
    a.dissolved = true;
}


void Mesh::assignQuadsToRegions(){
    quadRegions.init( quads.size() );
    for (Edge e: edges) {
        if (!e.isCut) quadRegions.fuse( qiOf(e), qjOf(e) );
    }
    std::cout << "* regions: " << nRegions() << std::endl;

}

void Mesh::populateNodes(){
    nodes.clear();
    for (Vert&v: verts) v.ni = -1;

    std::vector<int> archValency( verts.size(), 0);
    for (Edge e: edges) {
        if (e.isCut) {
            archValency[ viOf(e) ]++;
            archValency[ vjOf(e) ]++;
        }
    }
    for (int vi=0; vi<verts.size(); vi++ ){
        if ((archValency[vi] != 2 && archValency[vi] != 0)||(!verts[vi].isRegular())) {
            verts[ vi ].ni = nodes.size();
            nodes.push_back( Node{ vi, archValency[vi] } );
        }
    }
    std::cout << "* nodes: " << nodes.size() << std::endl;
}

void Mesh::createArchFrom(Half h){
    Arc arc;
    arc.ni = vert(h).ni;
    arc.nj = -1;
    arc.dissolved = false;
    arc.ei.clear();


    while (arc.nj==-1) {
        int vj = opposite(h).vi;
        arc.nj = verts[vj].ni; // -1 if not a node
        assert( edge(h).isCut );
        assert( !isMarked(edge(h)) );

        mark( edge(h) );
        arc.ei.push_back( h.ei );

        h = ahead( h );
    }

    arcs.push_back(arc);
    nodes[arc.ni].valency++;
    nodes[arc.nj].valency++;

}

bool Mesh::dissolveArches(){
    //std::cout<<"PASS (" << arcs.size() << " arcs left)" << std::endl;
    std::sort(arcs.begin(),arcs.end()); // prioritization
    for (Arc& a:arcs) {
        if (canBeDissolved(a)) dissolve(a);
    }
    return compressArcsVector();
}

bool Mesh::compressArcsVector(){
    std::vector<Arc> surviving;
    bool anyoneDied = false;
    for (Arc& a:arcs) {
        if (!a.dissolved) surviving.push_back(a);
        else anyoneDied = true;
    }
    if (!anyoneDied) return false;
    arcs = surviving;
    return true;
}

void Mesh::populateArches(){
    unmarkAll();
    arcs.clear();
    for (Node& n:nodes) n.valency = 0;

    for (Node n:nodes) {
        int hi = verts[ n.vi ].hi;
        int start = hi;
        do {
            hi = pivotCW(hi);
            Half h = halfs[hi];
            if (edges[h.ei].isCut && !isMarked(edge(h))) createArchFrom(h);
        } while (hi != start);
    }

    for (Arc& a:arcs) a.priority = evaluatePriority(a);
    std::cout << "* arcs: " << arcs.size() << std::endl;
}

void Mesh::print(const Arc& a) const{
    std::cout << "Arch from "<< nodes[a.ni].vi << " (" << nodes[a.ni].valency << ")";
    std::cout << " to " << nodes[a.nj].vi << " (" << nodes[a.nj].valency << "): ";
}

void Mesh::print(const Node& n) const{
    std::cout << "+node "<< n.vi << ": arch-val="<< n.valency << ",edge-val=" << verts[n.vi].valency << "\n";
}
void Mesh::printNodes() const{
    std::cout << "Nodes["<< nodes.size() << "]:\n";
    for (Node n:nodes) print(n);
}
