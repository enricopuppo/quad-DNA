#pragma once
#include<vector>
#include<string>

class PrimaryStructure;
class DNA;

struct Chunk{
    int nodeId;
    int nodeVal;
    int spanId;
};

class RNA{
public:

    std::string toString() const;
    bool fromString(const std::string&);

    int evaluateGenus() const;

    void markForDebug(); // this makes it easiers to follow <=> DNA, but breakes it

private:
    std::vector<Chunk> chunks;
    std::vector<int> spanLenghts;
    void clear();

    int edgeLenghtOf(const Chunk &c) const;
    int vertLenghtOf(const Chunk &c) const;

    int maxSpanId() const { return chunks.size()/2;}
    int otherIndexOfSameSpan(int i) const;
    // to preserve the mapping
    //std::vector<int> origVi;

friend PrimaryStructure translation(const RNA &);
friend RNA reverseTranslation(const PrimaryStructure &);
friend RNA transcription(const DNA &);  // the work of the ribosome :)
friend DNA reverseTranscription(const RNA &);

};

PrimaryStructure translation(const RNA &);  // the work of the ribosome :)
RNA reverseTranslation(const PrimaryStructure &);

