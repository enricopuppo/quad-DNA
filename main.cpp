#include <iostream>
#include <fstream>

#include "mesh.h"
#include "primary_structure.h"
#include "RNA.h"
#include "DNA.h"
#include "gl_viewer.h"

const char *downArrowAsciiArt = "     v v v \n"; // "||\n  \\/\n"; ↓

std::string filename;
std::ofstream logFile;

DNA encode(Mesh &m)
{
    logFile << "------- Encoding phase:\n";
    // RNA encode(Mesh& m){
    m.computeMotorcycleGraphCuts();
    // m.exportCutsAsObj(filename+".cut.MCG.obj");

    m.dissolveCutsToMakeOneDisk();
    m.exportCutsAsObj(filename + ".cut.final.obj");

    logFile << "Secondary structure: " << m.statsString() << std::endl;

    PrimaryStructure prim = reversePrefolding(m);
    prim.setLogFile(logFile);
    prim.exportAsObj(filename + ".uv.border.obj");
    logFile << downArrowAsciiArt << "Primary structure: " << prim.toString(3000) << std::endl;

    RNA rna = reverseTranslation(prim);
    // rna.markForDebug();
    logFile << downArrowAsciiArt << "RNA: " << rna.toString() << std::endl;

    DNA dna = reverseTranscription(rna);
    logFile << downArrowAsciiArt << dna.toStringArt("DNA: ") << std::endl;

    // to remap vertices on the mesh...
    auto x = translation(transcription(dna));
    // auto x = translation( rna );
    logFile << downArrowAsciiArt << "primary structure: " << x.toString() << std::endl;
    logFile << "-------- End of encoding phase --------\n\n";

    prim.renameVerticesAs(x);
    prim.prefolding(&m);
    m.renameVertices(prim.origToNewVi);

    return dna;
}

Mesh decode(const DNA &dna)
{
    logFile << "------- Decoding phase:\n";

    RNA rna = transcription(dna);
    logFile << downArrowAsciiArt << "RNA: " << rna.toString() << std::endl;

    PrimaryStructure prim = translation(rna);
    logFile << downArrowAsciiArt << "Primary structure: " << prim.toString() << std::endl;

    std::vector<Quad> quads = prim.prefolding();

    logFile << downArrowAsciiArt << "Secondary structure: " << quads.size() << " quads" << std::endl;

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

    logFile << downArrowAsciiArt << "Secondary structure: " << quads.size() << " quads" << std::endl;
    logFile << "-------- End of decoding phase --------\n\n";

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
    logFile.open(filename + ".log");
    mi.setLogFile(logFile);

    if (!mi.importObj(filename))
    {
        logFile << "ABORT - File not opened or mesh not closed\n\n\n\n\n";
        logFile.close();
        return 1;
    }

    auto dna = encode(mi);

    Mesh mo = decode(dna);
    Mesh uv = decodeUV(dna);

    mo.assignEmbeddingAs(mi);
    mo.setLogFile(logFile);
    mo.exportObj(filename + ".reconstructed.obj", uv);
    uv.setLogFile(logFile);
    uv.exportObj(filename + ".uvmap.obj");

    // GlViewer viewer;
    // viewer.openWindow();
    // viewer.renderMesh(mo);
    // viewer.waitForClose();

    logFile.close();
    return 0;
}
