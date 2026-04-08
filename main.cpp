#include <iostream>

#include "mesh.h"
#include "primary_structure.h"
#include "RNA.h"
#include "DNA.h"
#include "gl_viewer.h"

const char *downArrowAsciiArt = "     v v v \n"; // "||\n  \\/\n"; ↓

std::string filename;

DNA encode(Mesh &m)
{
    // RNA encode(Mesh& m){
    m.computeMotorcycleGraphCuts();
    // m.exportCutsAsObj(filename+".cut.MCG.obj");

    m.dissolveCutsToMakeOneDisk();
    m.exportCutsAsObj(filename + ".cut.final.obj");

    std::cout << "\n\nSecondary structure: " << m.statsString() << std::endl;

    PrimaryStructure prim = reversePrefolding(m);
    prim.exportAsObj(filename + ".uv.border.obj");
    std::cout << downArrowAsciiArt << "Primary structure: " << prim.toString(3000) << std::endl;

    RNA rna = reverseTranslation(prim);
    // rna.markForDebug();

    std::cout << downArrowAsciiArt << "encode RNA: " << rna.toString() << std::endl;

    DNA dna = reverseTranscription(rna);

    std::cout << downArrowAsciiArt << dna.toStringArt("DNA: ") << std::endl;

    // to remap vertices on the mesh...
    auto x = translation(transcription(dna));
    // auto x = translation( rna );
    std::cout << downArrowAsciiArt << "primary structure: " << x.toString() << std::endl;

    prim.renameVerticesAs(x);
    prim.prefolding(&m);
    m.renameVertices(prim.origToNewVi);

    return dna;
}

Mesh decode(const DNA &dna)
{

    RNA rna = transcription(dna);
    std::cout << downArrowAsciiArt << "decode RNA: " << rna.toString() << std::endl;

    PrimaryStructure prim = translation(rna);
    std::cout << downArrowAsciiArt << "Primary structure: " << prim.toString() << std::endl;

    std::vector<Quad> quads = prim.prefolding();

    std::cout << downArrowAsciiArt << "Secondary structure: " << quads.size() << " quads" << std::endl;
    std::cout << "\n\n";

    Mesh m;
    m.assignConnectivity(quads);
    // no embedding!

    return m;
}

Mesh decodeUV(const DNA &dna)
{

    RNA rna = transcription(dna.nonCanonicalVersion());

    PrimaryStructure prim = translation(rna);

    prim.splitAllCuts();
    std::vector<Quad> quads = prim.prefolding();

    std::cout << downArrowAsciiArt << "Secondary structure: " << quads.size() << " quads" << std::endl;
    std::cout << "\n\n";

    Mesh m;
    m.assignConnectivity(quads);
    m.assignEmbeddig(prim.uvMap);
    m.normalizeIn01();

    return m;
}

int main(int argc, char *argv[])
{

    assert(argc == 2);
    filename = argv[1];

    Mesh mi;

    mi.importObj(filename);

    auto dna = encode(mi);

    Mesh mo = decode(dna);
    Mesh uv = decodeUV(dna);

    mo.assignEmbeddingAs(mi);
    mo.exportObj(filename + ".reconstructed.obj", uv);
    uv.exportObj(filename + ".uvmap.obj");

    GlViewer viewer;
    viewer.openWindow();
    viewer.renderMesh(mo);
    viewer.waitForClose();

    return 0;
}
