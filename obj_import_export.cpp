/* implementation of all OBJ imporer/exporter methods */

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include "mesh.h"
#include "primary_structure.h"

class OutputObjFile
{
    std::ofstream file;
    std::ostream *logFile = &std::cout; // for diagnostics, not necessarily the same as file

public:
    void setLogFile(std::ostream &out) { logFile = &out; }
    std::ostream &getLogFile() { return *logFile; }
    const std::ostream &getLogFile() const { return *logFile; }

    OutputObjFile(const std::string &filename, std::ostream &out) : file(filename), logFile(&out)
    {
        (*logFile) << "OBJ file: " << filename << std::endl;
        if (!isOpen())
        {
            std::cerr << "* Failed to write on file" << std::endl;
            return;
        }
        file << "obj\n";
    }
    bool isOpen() const
    {
        return file.is_open();
    }
    void savePos(glm::vec3 pos)
    {
        file << "v " << pos.x << " " << pos.y << " " << pos.z << "\n";
    }
    void savePos(UvPos pos, float z = 0)
    {
        file << "v " << pos.u << " " << pos.v << " " << z << "\n";
    }
    void saveUv(glm::vec3 pos)
    {
        file << "vt " << pos.x << " " << pos.y << "\n";
    }
    void saveEdge(int a, int b, int vn)
    {
        file << "f " << (a + 1) << " " << (b + 1) << " " << (b + vn + 1) << " " << (a + vn + 1) << "\n";
    }
    void saveQuad(Quad q)
    {
        file << "f";
        for (int w = 0; w < 4; w++)
            file << " " << (q.vi[w] + 1);
        file << "\n";
    }
    void saveQuad(Quad qa, Quad qb)
    {
        file << "f";
        for (int w = 0; w < 4; w++)
            file << " " << (qa.vi[w] + 1) << "/" << (qb.vi[w] + 1);
        file << "\n";
    }

    ~OutputObjFile()
    {
        if (file.is_open())
            file.close();
    }
};

class InputObjFile
{
    std::ifstream file;
    std::ostream *logFile = &std::cout; // for diagnostics, not necessarily the same as file

public:
    void setLogFile(std::ostream &out) { logFile = &out; }
    std::ostream &getLogFile() { return *logFile; }
    const std::ostream &getLogFile() const { return *logFile; }
    InputObjFile(const std::string &filename, std::ostream &out) : file(filename), logFile(&out)
    {
        (*logFile) << "Loading OBJ file: " << filename << std::endl;
        if (!isOpen())
        {
            std::cerr << "* Failed to read file" << std::endl;
            return;
        }
        std::string line;
        std::getline(file, line);
        if (line == "obj")
        {
            std::cerr << "* Not an OBJ file?" << std::endl;
            file.close();
        }
    }
    bool isOpen() const
    {
        return file.is_open();
    }
    void readVertsAndQuads(std::vector<Vert> &verts, std::vector<Quad> &quads)
    {
        std::string line;

        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;

            std::istringstream ss(line);
            std::string prefix;
            ss >> prefix;

            if (prefix == "v")
            {
                Vert v;
                ss >> v.pos.x >> v.pos.y >> v.pos.z;
                verts.push_back(v);
            }
            else if (prefix == "f")
            {
                Quad q;
                for (int i = 0; i < 4; ++i)
                {
                    std::string token;
                    ss >> token;

                    if (token.empty())
                    {
                        std::cerr << "Invalid face (not a quad)" << std::endl;
                        break;
                    }

                    std::istringstream ts(token);
                    std::string indexStr;
                    std::getline(ts, indexStr, '/');

                    int index = std::stoi(indexStr);

                    q.vi[i] = index - 1; // OBJ indices are 1-based
                }
                quads.push_back(q);
            }
            // Ignore everything else (vt, vn, etc.)
        }
    }

    ~InputObjFile()
    {
        if (file.is_open())
            file.close();
    }
};

bool Mesh::importObj(const std::string &filename)
{
    verts.clear();
    quads.clear();

    InputObjFile file(filename, *logFile);
    if (!file.isOpen())
        return false;

    file.readVertsAndQuads(verts, quads);

    populateHalfs();
    if (!isClosed())
        return false; // std::cout << "WARNING mesh not closed!" << std::endl;
    (*logFile) << "* " << quads.size() << " quad faces and " << verts.size() << " verts." << std::endl;

    populateEdges();
    (*logFile) << "* " << edges.size() << " edges." << std::endl;

    updateVertexValencies();
    (*logFile) << "* " << countIrregulars() << " irregular vertices" << std::endl;

    return true;
}

bool Mesh::exportCutsAsObj(const std::string &filename) const
{
    (*logFile) << "Saving cuts in ";
    OutputObjFile file(filename, *logFile);
    if (!file.isOpen())
        return false;

    const glm::vec3 smallDelta(0.001, 0.00, 0.001);
    for (Vert v : verts)
        file.savePos(v.pos);
    for (Vert v : verts)
        file.savePos(v.pos + smallDelta);

    for (Edge e : edges)
    {
        if (e.isCut)
            file.saveEdge(viOf(e), vjOf(e), verts.size());
    }
    return true;
}

bool Mesh::exportObj(const std::string &filename) const
{
    (*logFile) << "Saving mesh in ";
    OutputObjFile file(filename, *logFile);
    if (!file.isOpen())
        return false;

    for (Vert v : verts)
        file.savePos(v.pos);
    for (Quad q : quads)
        file.saveQuad(q);

    return true;
}

bool Mesh::exportObj(const std::string &filename, const Mesh &uvMap) const
{
    (*logFile) << "Saving mesh in ";
    OutputObjFile file(filename, *logFile);
    if (!file.isOpen())
        return false;

    int nq = quads.size();
    assert(nq == uvMap.quads.size());

    for (Vert v : verts)
        file.savePos(v.pos);
    for (Vert v : uvMap.verts)
        file.saveUv(v.pos);
    for (int qi = 0; qi < nq; qi++)
        file.saveQuad(quads[qi], uvMap.quads[qi]);

    return true;
}

bool PrimaryStructure::exportAsObj(const std::string &filename) const
{
    (*logFile) << "Saving 2D boundary in ";
    OutputObjFile file(filename, *logFile);
    if (!file.isOpen())
        return false;

    for (int i = 0; i < sequence.size(); i++)
    {
        file.savePos(sequence[i].pos(), 2.0 * i / sequence.size());
    }
    for (int i = 0; i < sequence.size(); i++)
    {
        file.savePos(sequence[i].pos(), 2.0 * i / sequence.size() + 0.01);
    }

    int i = sequenceStart;
    int prev = -1;
    do
    {
        AminoAcid a = sequence[i];
        if (prev != -1)
            file.saveEdge(prev, i, sequence.size());
        prev = i;
        i = a.next;
    } while (i != sequenceStart);
    file.saveEdge(prev, sequenceStart, sequence.size());

    return true;
}
