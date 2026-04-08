#pragma once
#include <vector>

class UnionFind{
public:
    void init(int nRegions);
    void fuse(int regionA, int regionB);
    int areSameCluster(int regionA, int regionB);
    int nClusters() const { return nDistinct; };

private:
    int find(int i);
    std::vector<int> boss;
    int nDistinct;
};
