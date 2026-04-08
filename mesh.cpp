#include <map>
#include <glm/common.hpp>

#include "mesh.h"
#include "uv.h"


BoundingCube Mesh::computeBoundingCube() const{
    if (verts.empty()) return BoundingCube();
    glm::vec3 bbMin, bbMax;
    bbMin = bbMax = verts[0].pos;
    for (const Vert &v:verts) {
        bbMin = glm::min( bbMin, v.pos );
        bbMax = glm::max( bbMax, v.pos );
    }
    BoundingCube res;
    res.center = (bbMin+bbMax)/2.0f;
    res.radius = glm::max(bbMax.x-bbMin.x, glm::max( bbMax.y-bbMin.y, bbMax.z-bbMin.z )) / 2.0f;

    return res;

}

void Mesh::normalizeIn01(){
    if (verts.empty()) return;
    glm::vec3 bbMin, bbMax;
    bbMin = bbMax = verts[0].pos;
    for (const Vert &v:verts) {
        bbMin = glm::min( bbMin, v.pos );
        bbMax = glm::max( bbMax, v.pos );
    }
    float dim = glm::max( bbMax.x-bbMin.x, bbMax.y-bbMin.y);
    for (Vert &v:verts) {
        v.pos = (v.pos - bbMin)/dim;
    }
}
void Mesh::populateVerts(){
    verts.clear();
    int maxV = 0;
    for (Quad q:quads) {
        for (int vi:q.vi)
        maxV = std::max( maxV, vi);
    }
    verts.resize(maxV);
}

std::string Mesh::statsString() const{
    std::string res;
    int cutEdges = 0;
    for (const Edge& e:edges) if (e.isCut) cutEdges++;
    res += std::to_string(quads.size()) + " quads ("
           + std::to_string(verts.size()) + " distinct verts,  "
           + std::to_string(cutEdges) + " cut edges)";
    return res;
}

/*static void print(std::vector<int> v, const char *title){
    std::cout << title << ": (";
    for (int a:v) std::cout << a << " ";
    std::cout << ")" << std::endl;
}*/

void Mesh::renameVertices(const std::vector<int> newVi){

    //::print(newVi,"RENAMING");
    for (Half &h:halfs) h.vi = newVi[ h.vi ];
    for (Quad &q:quads) for (int w=0; w<4; w++)  q.vi[w] = newVi[ q.vi[w] ];
    for (Node &n:nodes) n.vi = newVi[ n.vi ];

    auto copy = verts;
    // TODO: remove checks!
    for (int i=0; i<newVi.size(); i++) if (newVi[i]>=0 && newVi[i]<verts.size()) verts[ newVi[i] ] = copy[ i ];
}

typedef std::pair<int,int> ViPair;
static ViPair ordered(ViPair v) { if (v.first>v.second) return v; else return ViPair(v.second,v.first); }

void Mesh::populateHalfs(){

    halfs.clear();

    std::map<ViPair,int> map;

    for (int qi = 0; qi<quads.size(); qi++) {
        Quad &q = quads[qi];
        int hi[4];

        int nh = halfs.size();
        hi[0] = nh+0;
        hi[1] = nh+1;
        hi[2] = nh+2;
        hi[3] = nh+3;
        halfs.resize( nh+4 );

        for (int i=0; i<4; i++) {
            int vi = q.vi[i];
            int vj = q.vi[(i+1)%4];

            Half &h = halfs[ hi[i] ];
            h.vi = vi;
            h.qi = qi;
            h.next = hi[ (i+1)%4 ];

            quads[qi].hi = verts[vi].hi = hi[i];

            ViPair pji(vj,vi);

            auto found = map.find(pji);
            if (found==map.end()) {
                h.opp = -1;
                ViPair pij(vi,vj);
                map[pij] = hi[i];
            } else{
                h.opp = found->second;
                halfs[ h.opp ].opp = hi[i];
            }

        }
    }
}

void Mesh::clearCuts(){
    for (Edge& e:edges) e.isCut = false;
}

void Mesh::populateEdges(){
    edges.clear();
    unmarkAll();

    for (int hi=0; hi<halfs.size(); hi++) {
        Half &h = halfs[hi];
        if (!isMarked(h)) {
            Edge e;
            e.hi = hi;
            h.ei = edges.size();
            edges.push_back( e );
            mark( h );
            mark( opposite(h) );
        } else {
            h.ei = opposite( h ).ei;
        }
    }

}

bool Mesh::isClosed() const{
    for (const Half &h : halfs){
        if (h.opp == -1) return false;
    }
    return true;
}

void Mesh::updateVertexValencies(){
    for (Vert &v: verts) {
        v.valency = 0;
    }
    /*for (Edge e: edges) {
        verts[ viOf(e) ].valency ++;
        verts[ vjOf(e) ].valency ++;
    }*/
    for (Quad q: quads) {
        verts[ q.vi[0] ].valency ++;
        verts[ q.vi[1] ].valency ++;
        verts[ q.vi[2] ].valency ++;
        verts[ q.vi[3] ].valency ++;
    }
}

int Mesh::countIrregulars() const{
    int res = 0;
    for (Vert v: verts) {
        if (!v.isRegular()) res++;
    }
    return res;
}

void Mesh::assignEmbeddingAs(const Mesh& other){
    verts.resize(other.verts.size());
    for (int i=0;  i<verts.size(); i++) {
        verts[i].pos = other.verts[i].pos;
    }
}

void Mesh::assignEmbeddig(const std::vector<UvPos>& uvmap){
    verts.resize(uvmap.size());

    updateVertexValencies();
    // z as creation time
    std::vector<float> zQuad(quads.size());
    for (int i=0; i<quads.size(); i++) zQuad[i] = 2.0*float(i)/quads.size();

    // smooth z
    std::vector<float> zVert(verts.size());
    scatterVertAttr(zQuad,zVert);
    gatherQuadAttr(zVert,zQuad);
    scatterVertAttr(zQuad,zVert);

    for (int i=0;  i<verts.size(); i++) {
        verts[i].pos = glm::vec3( uvmap[i].u, uvmap[i].v, zVert[i]) ;
    }
}


int Mesh::findStartingCutHalf() const{
    /* first attempt: a deadend ("baffo") */
    for (int hi=0; hi<halfs.size(); hi++) {
        Half h = halfs[hi];
        if (edge(h).isCut) {
            int ni = vertEnd(h).ni;
            if (ni>=0) {
                int av = nodes[vertEnd(h).ni].valency;
                if (av == 1)  return hi;
            }
        }
    }

    /* 2st attempt: a irregular point*/
    for (int hi=0; hi<halfs.size(); hi++) {
        Half h = halfs[hi];
        if (edge(h).isCut) {
            if (!vertEnd(h).isRegular()) return hi;
        }
    }

    // no irregular vert?
    int hi = halfs[verts[0].hi].opp; // half pointing to vert 0
    assert( edge(halfs[hi]).isCut );
    return hi;
}

int Mesh::followCut(int hi, int & valency) const{
    Half h = halfs[hi];
    assert(edge(h).isCut);
    valency = 1;
    hi = next(hi);
    while (!edge(halfs[hi]).isCut){
        hi = pivotCW(hi);
        valency++;
    }
    return hi;
}

void Mesh::scatterVertAttr(const std::vector<float>& perQuad, std::vector<float>& perVert) const{
    assert(perQuad.size()==quads.size());
    assert(perVert.size()==verts.size());
    for (float& pv:perVert) pv = 0;
    int k=0;
    for (Quad q:quads) {
        for (int w=0;w<4;w++) perVert[q.vi[w]] += perQuad[k];
        k++;
    }
    k=0;
    for (const Vert& v:verts) {
        perVert[k++] /= v.valency;
    }
}

void Mesh::gatherQuadAttr(const  std::vector<float>& perVert, std::vector<float>& perQuad ) const{
    assert(perQuad.size()==quads.size());
    assert(perVert.size()==verts.size());
    int k=0;
    for (Quad q:quads) {
        float sum=0;
        for (int w=0;w<4;w++) sum+= perVert[q.vi[w]];
        perQuad[k++] = sum/4;
    }
};


