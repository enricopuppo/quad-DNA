#include "union_find.h"

void UnionFind::init(int n){
    nDistinct = n;
    boss.resize(n);
    for (int i=0; i<n; i++) boss[i] = i;
}

static bool arbitraryChoice(int i, int j){
    return (i+j)%2 == 0;
}
void UnionFind::fuse(int i, int j){
    int bi = find(i);
    int bj = find(j);
    if (bi==bj) return;
    nDistinct--;
    if (arbitraryChoice(i,j)) boss[bi] = bj;
    else boss[bj] = bi;
}

int UnionFind::areSameCluster(int i, int j){
    return (find(i)==find(j));
}

int UnionFind::find(int i){
    while (1) {
        int j = boss[i];
        if (j == i) return i;
        boss[i] = j; // path shortening
        i = j;
    }
}
