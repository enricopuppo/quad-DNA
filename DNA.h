#pragma once
#include <string>
#include <vector>

class RNA;
class Chunk;

struct BasePair{
    // in this precise order:
    int base[2]; // the two turns
    std::vector<int> portsIn;
    bool closedRound = false;
    bool portOut = false;
    int openRounds = 0;
    int len;
    // methods
    void setZero();
    bool isCanon() const;
    BasePair splitCanonPad();
};

class DNA
{
public:
    std::string toString() const;
    std::string toStringFat() const;
    std::string toStringArt(const char* prefix="") const;
    bool fromString(const std::string&);

    DNA canonicalVersion() const;
    DNA nonCanonicalVersion() const;

private:

    std::vector<BasePair> helix; // the double helix. cyclic!

    const BasePair& basePair(int index) const { return helix[index % helix.size()]; /* it's cyclcic */ }
    BasePair& basePair(int index)  { return helix[index % helix.size()];}

    struct Site{
        int index;  // where in the double helix
        int side;   // which side of the double helix. Determines the verse: 0: increase. 1: decrease
        bool endReached = false;
        bool operator != (const Site& other) const{ return index!=other.index || side!=other.side; }
        std::string toString() const;
        int indexOfVi;
        bool isOver() const { return side==1 && index==0; /*endReached;*/ }
    };

    Site step(Site s) const;
    Site stepForward(Site s) const;
    Site stepBackward(Site s) const;

    bool findMatchingOpenParOf(int &helixPos) const; // true if it ends in another open

    Chunk chunkAt(Site s) const;

    friend RNA transcription(const DNA &);  // the work of the ribosome :)
    friend DNA reverseTranscription(const RNA &);

    // teleporting locations
    Site exitOfInPort(Site i, int which) const;
    Site exitOfLastInPort(Site i) const;
    Site exitOfOutPort(Site i) const;

    std::pair<int,int> ithInPort(int i) const;
    int ithOutPort(int i) const;
    int rankOfOutPort(int helixPos) const;

};

RNA transcription(const DNA &);
DNA reverseTranscription(const RNA &);
