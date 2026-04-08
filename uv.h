#pragma once

enum UvDir{
    N,E,S,W
};

inline UvDir operator + (UvDir a, int valency){
    return UvDir((int(a)+valency+2)%4);
}

struct UvPos{
    int u, v;
    UvPos(int u,int v):u(u),v(v){}
    UvPos():UvPos(0,0){}
    UvPos operator + (UvDir b){
        switch(b) {
        default:
        case N: return UvPos( u ,v+1);
        case S: return UvPos( u ,v-1);
        case E: return UvPos(u+1, v );
        case W: return UvPos(u-1, v );
        }
    }
    UvPos operator - (UvDir b){
        switch(b) {
        default:
        case N: return UvPos( u ,v-1);
        case S: return UvPos( u ,v+1);
        case E: return UvPos(u-1, v );
        case W: return UvPos(u+1, v );
        }
    }

    bool operator == (const UvPos& b) const { return u==b.u && v == b.v; }
};
