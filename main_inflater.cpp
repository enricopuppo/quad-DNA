#include <iostream>
#include <fstream>

#include "mesh.h"
#include "gl_viewer.h"
#include <math.h>

const char *downArrowAsciiArt = "     v v v \n"; // "||\n  \\/\n"; ↓

std::string filename;
std::ofstream logFile;

Mesh inflate(const Mesh &m)
{
    const float angle_THRESHOLD = M_PI/6.0; // in radians, to determine which irregular vertices to constrain during inflation
    const float K_diag = 1.0 / sqrt(2); // stiffness of diagonal edges, relative to regular edges
    Mesh inflated = m;
    // set initial constraints as irregular vertices with a sufficiently large angle defect/excess
    std::vector<int> irregulars = m.getIrregulars();
    std::vector<int> constraints;
    for (int vi : irregulars)
        if (fabs(m.angleDefect(vi)) > angle_THRESHOLD)
            constraints.push_back(vi);



    return inflated;
}

int main(int argc, char *argv[])
{

    assert(argc == 2);
    filename = argv[1];

    Mesh mi;
    // logFile.open(filename + ".log");
    // mi.setLogFile(logFile);

    if (!mi.importObj(filename))
    {
        std::cout << "ABORT - File not opened or mesh not closed\n";
        // logFile << "ABORT - File not opened or mesh not closed\n";
        // logFile.close();
        return 1;
    }

    Mesh mo = inflate(mi);


    mo.exportObj(filename + ".inflated.obj");

    GlViewer viewer;
    viewer.openWindow();
    viewer.renderMesh(mo);
    viewer.waitForClose();

    logFile.close();
    return 0;
}
