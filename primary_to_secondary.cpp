#include <iostream>
#include <map>
#include "primary_structure.h"
#include "mesh.h"

// should be fields...
static std::map<Tortoise, int> redFlags;
static const Mesh *refMesh = NULL;

bool isRedFlagged(const Tortoise &a)
{
    auto found = redFlags.find(a);
    if (found == redFlags.end())
        return false;
    return found->second > 0;
}

void addRedFlag(const Tortoise &a)
{
    auto found = redFlags.find(a);
    if (found == redFlags.end())
        redFlags[a] = 1;
    else
        redFlags[a] = found->second + 1;
}

void removeRedFlag(const Tortoise &a)
{
    redFlags[a] = redFlags[a] - 1;
    // remove from map if zero?
}

bool PrimaryStructure::isRedFlagged(const AminoAcid &a) const
{
    // navigate to opposite corner
    Tortoise t = a.tortoise;
    t.turn(1);
    t.stepForward();
    t.turn(1);
    t.stepForward();
    t.turn(1);
    return ::isRedFlagged(t);
}

void PrimaryStructure::addRedFlags(const AminoAcid &a)
{
    if (a.isConcave())
    {
        Tortoise t = a.tortoise;
        for (int k = 2; k <= a.val - 1; k++)
        {
            addRedFlag(t);
            t.turn(3);
        }
    }
}

void PrimaryStructure::removeRedFlags(const AminoAcid &a)
{
    if (a.isConcave())
    {
        Tortoise t = a.tortoise;
        for (int k = 2; k <= a.val - 1; k++)
        {
            removeRedFlag(t);
            t.turn(3);
        }
    }
}

void PrimaryStructure::decreaseValency(AminoAcid &a)
{
    a.val--;
    assert(a.val > 0);
}

void PrimaryStructure::assignValency(AminoAcid &a, int newVal)
{
    a.val = newVal;
}

void PrimaryStructure::removePairFromSequence(AminoAcid &a, AminoAcid &b)
{
    AminoAcid &before = sequence[a.prev];
    AminoAcid &after = sequence[b.next];

    int toA = before.next;
    int toB = after.prev;
    if (sequenceStart == toA || sequenceStart == toB)
        sequenceStart = b.next;

    before.next = b.next;
    after.prev = a.prev;
}

Quad PrimaryStructure::bonkRetract(AminoAcid &a, AminoAcid &b, AminoAcid &c, AminoAcid &d)
{
    removeRedFlags(a);
    removeRedFlags(d);

    // std::cout << "BONK!" << std::endl;
    Quad q;
    q.vi[0] = a.vi;
    q.vi[1] = b.vi;
    q.vi[2] = c.vi;
    q.vi[3] = d.vi;

    removePairFromSequence(b, c);

    decreaseValency(a);
    decreaseValency(d);
    d.tortoise.turn(3);

    addRedFlags(a);
    addRedFlags(d);

    if (refMesh)
    {
        a.hiOnMesh = refMesh->opposite(refMesh->next(c.hiOnMesh));
    }
    return q;
}

Quad PrimaryStructure::cornerRetract(AminoAcid &a, AminoAcid &b, AminoAcid &c)
{
    // std::cout << "CORNER!" << std::endl;
    removeRedFlags(a);
    removeRedFlags(c);

    b.tortoise.stepBackward();
    b.tortoise.turn(1);
    b.tortoise.stepForward();

    c.tortoise.turn(3);

    Quad q;
    int vi = newVertIndex(b.pos());
    q.vi[0] = a.vi;
    q.vi[1] = b.vi;
    q.vi[2] = c.vi;
    q.vi[3] = vi;
    b.vi = vi;
    assignValency(b, 3);
    decreaseValency(a);
    decreaseValency(c);

    addRedFlags(a);
    addRedFlags(b);
    addRedFlags(c);

    if (refMesh)
    {
        a.hiOnMesh = refMesh->opposite(refMesh->next(refMesh->next(b.hiOnMesh)));
        b.hiOnMesh = refMesh->opposite(refMesh->next(b.hiOnMesh));
        int viOnMesh = refMesh->halfs[b.hiOnMesh].vi;
        origToNewVi[viOnMesh] = vi;
        // if (vi==1000) std::cout<<"FOUND "<< viOnMesh << " : "<< vi << std::endl;
    }

    return q;
}

Quad PrimaryStructure::finalQuad(AminoAcid &a, AminoAcid &b, AminoAcid &c, AminoAcid &d) const
{
    // std::cout << "FINAL!" << std::endl;
    Quad q;
    q.vi[0] = a.vi;
    q.vi[1] = b.vi;
    q.vi[2] = c.vi;
    q.vi[3] = d.vi;
    return q;
}

void PrimaryStructure::setInitialRedFlags()
{
    redFlags.clear();
    int i = sequenceStart;
    do
    {
        const AminoAcid &a = sequence[i];
        addRedFlags(a);
        i = a.next;
    } while (i != sequenceStart);
}

std::vector<Quad> PrimaryStructure::prefolding(const Mesh *_refMesh)
{

    refMesh = _refMesh;
    if (refMesh == NULL)
        origToNewVi.clear();

    // std::cout << "START SHRINK " << toString() << "!!!\n";

    std::vector<Quad> quads;
    setInitialRedFlags();

    int i = sequenceStart;

    while (1)
    {
        AminoAcid &a = sequence[i];
        AminoAcid &b = sequence[a.next];
        AminoAcid &c = sequence[b.next];

        if (b.isConvex())
        {
            if (c.isConvex())
            {
                AminoAcid &d = sequence[c.next];
                if (d.isConvex())
                {
                    quads.push_back(finalQuad(a, b, c, d));
                    break;
                }
                else
                {
                    quads.push_back(bonkRetract(a, b, c, d));
                    i = a.prev;
                }
            }
            else
            {
                if (!isRedFlagged(b) && !a.isConvex())
                    quads.push_back(cornerRetract(a, b, c));
            }
        }
        // todo: DETECT DEADLOCKS (too many moves)
        i = a.next;
    }

    // std::cout << "2nd to 3rd: " << quads.size() << " quads and " << nv << " vertices procuded" << std::endl;
    // std::cout << "END SHRINK " << toString() << "!!!\n";
    return quads;
}
