#include "mesh.h"
#include <iostream>

// typedef int Motor;
struct Motor
{
    int hi;
    int patience = 0; // how many crashes before giving up. 0 for standard MCG
    Motor(int hi) : hi(hi) {}
    bool isAlive() const { return hi >= 0; }
    void kill() { hi = -1; }
    void crash()
    {
        if (patience-- <= 0)
            kill();
    }

    static Motor dead;
};

Motor Motor::dead = Motor(-1);

class MCG
{
public:
    Mesh &mesh;
    MCG(Mesh &m) : mesh(m) {}

    std::vector<Motor> motors; // each is one halfedge
    int turn = 0;

    void spawnAround(Vert &v)
    {
        mesh.mark(v);

        int hStart = v.hi;
        int hi = hStart;
        do
        {
            motors.push_back(hi);
            hi = mesh.pivotCW(hi);
        } while (hi != hStart);
    }

    void initMotors()
    {
        mesh.unmarkAll();
        for (Vert v : mesh.verts)
        {
            if (!v.isRegular())
                spawnAround(v);
        }

        if (motors.empty())
        {
            // no irreuglar? spawn around an arbitrary vert
            if (!mesh.verts.empty())
                spawnAround(mesh.verts[0]);
        }
        // std::cout << "MCG: spawned " << motors.size() << " motors!" << std::endl;
    }

    Motor &nextAliveMotor()
    {
        int initTurn = turn;
        do
        {
            turn = (turn + 1) % motors.size();
            if (motors[turn].isAlive())
                return motors[turn];
        } while (initTurn != turn);

        return Motor::dead;
    }

    bool oneStep()
    {
        Motor &mi = nextAliveMotor();
        if (!mi.isAlive())
            return false;
        oneStepMotor(mi);
        return true;
    }

    void oneStepMotor(Motor &m)
    {

        int ei = mesh.halfs[m.hi].ei;

        if (mesh.edges[ei].isCut)
        {
            m.kill(); // head on crush: sudden death
            return;
        }

        mesh.edges[ei].isCut = true;

        m.hi = mesh.ahead(m.hi);
        int vi = mesh.halfs[m.hi].vi;
        if (mesh.isMarked(mesh.verts[vi]))
        {
            m.crash();
        }
        mesh.mark(mesh.verts[vi]);
    }

    void exec()
    {
        initMotors();
        mesh.clearCuts();
        mesh.unmarkAll();
        turn = 0;
        while (oneStep())
            ;
    }
};

void Mesh::computeMotorcycleGraphCuts()
{
    MCG(*this).exec();
}
