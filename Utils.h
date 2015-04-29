/*
 * Utils.h
 *
 *  Created on: Jan 22, 2013
 *      Author: stef
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/graph_components.h>

using namespace lemon;
using namespace std;

#define DEBUG_PRINT(PREFIXSTR, EXPR) std::cout<<PREFIXSTR<<":\t"<<EXPR<<std::endl

#ifdef DEBUG_DATASET
#define DATASET_DEBUG_PRINT(EXPR) DEBUG_PRINT("[DATASET]", EXPR)
#else
#define DATASET_DEBUG_PRINT(EXPR) do {} while (0)
#endif

#ifdef DEBUG_TMP
#define TMP_DEBUG_PRINT(EXPR) DEBUG_PRINT("[TEMP]", EXPR)
#define TMP_DEBUG_PRINT_PAIR(x,y) cout<<x<<"\t"<<y<<endl
#define TMP_DEBUG_PRINT_TRIPLE(x,y,z) cout<<x<<"\t"<<y<<"\t"<<z<<endl
#define TMP_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define TMP_DEBUG_PRINT(EXPR) do {} while (0)
#define TMP_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define TMP_DEBUG_PRINT_TRIPLE(x,y,z) do {} while (0)
#define TMP_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_BP
#define BP_DEBUG_PRINT(EXPR) DEBUG_PRINT("[BASEPAIR]", EXPR)
#else
#define BP_DEBUG_PRINT(EXPR) do {} while (0)
#endif

#ifdef DEBUG_CLUSTER
#define CLUSTER_DEBUG_PRINT(EXPR) DEBUG_PRINT("[CLUSTER]", EXPR)
#else
#define CLUSTER_DEBUG_PRINT(EXPR) do {} while (0)
#endif

#ifdef DEBUG_FLOWGRAM
#define FLOWGRAM_DEBUG_PRINT(EXPR) DEBUG_PRINT("[FLOWGRAM]", EXPR)
#else
#define FLOWGRAM_DEBUG_PRINT(EXPR) do {} while (0)
#endif

#ifdef VERBOSE
//#define VERBOSE_PRINT(x) cout<<x<<endl
#define VERBOSE_PRINT(EXPR) DEBUG_PRINT("",EXPR)
#define VERBOSE_PRINT_LEVEL1(EXPR,LEVEL) if(1 <= LEVEL) DEBUG_PRINT("",EXPR)
#define VERBOSE_PRINT_LEVEL2(EXPR,LEVEL) if(2 <= LEVEL) DEBUG_PRINT("",EXPR)
//#define VERBOSE_PRINT_PAIR(x,y) cout<<x<<y<<endl
//#define VERBOSE_PRINT_TRIPLE(x,y,z) cout<<x<<y<<z<<endl
#define	VERBOSE_PRINT_PATH(path) printPath(path)
#else
#define VERBOSE_PRINT(EXPR) do {} while (0)
#define VERBOSE_PRINT_LEVEL1(EXPR,LEVEL) do {} while(0)
#define VERBOSE_PRINT_LEVEL2(EXPR,LEVEL) do {} while(0)
//#define VERBOSE_PRINT(x) do {} while (0)
//#define VERBOSE_PRINT_PAIR(x,y) do {} while (0)
//#define	VERBOSE_PRINT_TRIPLE(x,y, z) do {} while (0)
#define	VERBOSE_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_HOMOPOLYMER
#define	HOMOPOLYMER_PRINT_PATH(path) printPath(path)
#define HOMOPOLYMER_PRINT(EXPR) DEBUG_PRINT("[HOMOPOLYMER]", EXPR)
//#define HOMOPOLYMER_PRINT(x) cout<<x<<endl
//#define HOMOPOLYMER_PRINT_PAIR(x,y) cout<<x<<y<<endl
//#define HOMOPOLYMER_PRINT_TRIPLE(x,y,z) cout<<x<<y<<z<<endl
#else
#define	HOMOPOLYMER_PRINT_PATH(path) do {} while(0)
#define HOMOPOLYMER_PRINT(x) do {} while(0)
//#define HOMOPOLYMER_PRINT_PAIR(x,y)  do {} while(0)
//#define HOMOPOLYMER_PRINT_TRIPLE(x,y,z)  do {} while(0)
#endif

#ifdef DEBUG_ANTIBODY
#define ANTIBODY_DEBUG_PRINT(EXPR) DEBUG_PRINT("[ANTIBODY]", EXPR)
//#define ANTIBODY_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
//#define ANTIBODY_DEBUG_PRINT_TRIPLE(x,y,z) cout<<"FIXTIP:\t"<<x<<y<<z<<endl
#define ANTIBODY_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define ANTIBODY_DEBUG_PRINT(x) do {} while (0)
//#define ANTIBODY_DEBUG_PRINT_PAIR(x,y) do {} while (0)
//#define	ANTIBODY_DEBUG_PRINT_TRIPLE(x,y, z) do {} while (0)
#define ANTIBODY_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_COLOR
#define COLOR_DEBUG_PRINT(EXPR) DEBUG_PRINT("[COLOR]", EXPR)
//#define COLOR_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
//#define COLOR_DEBUG_PRINT_TRIPLE(x,y,z) cout<<"FIXTIP:\t"<<x<<y<<z<<endl
#define COLOR_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define COLOR_DEBUG_PRINT(x) do {} while (0)
//#define COLOR_DEBUG_PRINT_PAIR(x,y) do {} while (0)
//#define	COLOR_DEBUG_PRINT_TRIPLE(x,y, z) do {} while (0)
#define COLOR_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_CLIPTIP
#define CLIPTIP_DEBUG_PRINT(x) cout<<x<<endl
#define CLIPTIP_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#define CLIPTIP_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define CLIPTIP_DEBUG_PRINT(x) do {} while (0)
#define CLIPTIP_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define CLIPTIP_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_FIXTIP
#define FIXTIP_DEBUG_PRINT(x) cout<<"FIXTIP:\t"<<x<<endl
#define FIXTIP_DEBUG_PRINT_PAIR(x,y) cout<<"FIXTIP:\t"<<x<<y<<endl
#define FIXTIP_DEBUG_PRINT_TRIPLE(x,y,z) cout<<"FIXTIP:\t"<<x<<y<<z<<endl
#define FIXTIP_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define FIXTIP_DEBUG_PRINT(x) do {} while (0)
#define FIXTIP_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define FIXTIP_DEBUG_PRINT_TRIPLE(x,y, z) do {} while (0)
#define FIXTIP_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_BUBBLE
#define BUBBLE_DEBUG_PRINT(x) cout<<x<<endl
#define BUBBLE_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#define BUBBLE_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define BUBBLE_DEBUG_PRINT(x) do {} while (0)
#define BUBBLE_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define BUBBLE_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_DISTRIBUTE
#define DISTRIBUTE_DEBUG_PRINT(x) cout<<x<<endl
#define DISTRIBUTE_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#define DISTRIBUTE_DEBUG_PRINT_PATH(path) printPath(path)
#else
#define DISTRIBUTE_DEBUG_PRINT(x) do {} while (0)
#define DISTRIBUTE_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define DISTRIBUTE_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_ERRORCOR
//#define ERRORCOR_DEBUG_PRINT(x) cout<<x<<endl
#define ERRORCOR_DEBUG_PRINT(EXPR) DEBUG_PRINT("[ERRORCOR]", EXPR)
#define ERRORCOR_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#define ERRORCOR_DEBUG_PRINT_TRIPLE(x,y,z) cout<<x<<y<<z<<endl
#define ERRORCOR_DEBUG_PRINT_PATH(path) printPath(path)
#else
//#define ERRORCOR_DEBUG_PRINT(x) do {} while (0)
#define ERRORCOR_DEBUG_PRINT(EXPR) do {} while (0)
#define ERRORCOR_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define ERRORCOR_DEBUG_PRINT_TRIPLE(x,y,z) do {} while (0)
#define ERRORCOR_DEBUG_PRINT_PATH(path) do {} while (0)
#endif

#ifdef DEBUG_READPATH
#define READPATH_DEBUG_PRINT(x) cout<<x<<endl
#define READPATH_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#define READPATH_DEBUG_PRINT_PATH(path) printPath(path)
#define READPATH_DEBUG_PRINT_TRIPLE(x,y,z) cout<<"READPATH:\t"<<x<<y<<z<<endl
#else
#define READPATH_DEBUG_PRINT(x) do {} while (0)
#define READPATH_DEBUG_PRINT_PAIR(x,y) do {} while (0)
#define READPATH_DEBUG_PRINT_PATH(path) do {} while (0)
#define READPATH_DEBUG_PRINT_TRIPLE(x,y,z) do {} while (0)
#endif

#ifdef DEBUG_MERGE
#define MERGE_DEBUG_PRINT(x) cout<<x<<endl
#else
//#define DEBUG_PRINT(x) do {} while (0)
#endif

#ifdef DEBUG_MERGE
#define MERGE_DEBUG_PRINT_PAIR(x,y) cout<<x<<y<<endl
#else
#define DEBUG_PRINT_PAIR(x,y) do {} while (0)
#endif


#ifdef DEBUG_BUBBLE
#define _DEBUG_BUBBLE_PRINT_PATHS(src, mainStr, otherStr, sizeDiff, mainPath, lightPath) do {\
		cout<<"src: "<<g.id(src)<<endl;\
		cout<<"main str: "<<mainStr<<endl;\
		cout<<"other str:"<<otherStr<<endl;\
		cout<<"diff: "<<sizeDiff<<endl;\
		cout<<"main arc:"<<endl;\
		this->printPath(mainPath);\
		cout<<"light arc:"<<endl;\
		this->printPath(lightPath);\
		cout<<"-----------------------"<<endl;\
} while(0)
#else
#define DEBUG_BUBBLE_PRINT_PATHS(src, mainStr, otherStr, sizeDiff, mainPath, lightPath) do {} while (0)
#endif

#define VECTOR_PRINT(v) for(int iii = 0; iii < (int)v.size(); iii++) { cout<<v[iii]<<","; } cout<<endl;

inline string getRandomString(int len) {
	static const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	string str(len, 'x');
	for(int i = 0; i < len; i++) {
		str[i] = alphabet[rand() % sizeof(alphabet)-1];
	}
	return str;
}

inline string getSetAsString(set<int> s) {
	string str = "";
	for(set<int>::iterator sitr = s.begin(); sitr != s.end(); sitr++) {
		ostringstream ss;
		ss<< (*sitr);
		str += (ss.str()+",");
	}
	return str;
}

inline string getPairAsString(pair<int,int> p) {
	ostringstream ss2; ss2<< (p.first);
	ostringstream ss3; ss3<< (p.second);
	return "<" + ss2.str() + ", " + ss3.str() + ">";
}

inline string getVectorAsString(vector<int> v) {
	string str = "";
	for(int i = 0; i < (int)v.size(); i++) {
		str += (v[i] + ",");
	}
	return str;
}

inline string getMapAsString(map<int, pair<int,int> > m) {
	string str = "";
	for(map<int, pair<int,int> >::iterator sitr = m.begin(); sitr != m.end(); sitr++) {
		ostringstream ss1; ss1<< (sitr->first);
		str += ("[" + ss1.str() + ": " + getPairAsString(sitr->second) + "]" );
	}
	return str;
}
/**
 * return set that is intersection of two sets
 */
inline set<int> intersectionTwoSets(set<int> a, set<int> b) {
	set<int> tmp;
	for(set<int>::iterator aitr = a.begin(); aitr != a.end(); aitr++) {
		if(b.find(*aitr) != b.end()) tmp.insert(*aitr);
	}
	return tmp;
}

inline bool isIntersectionTwoSetsEmpty(set<int> a, set<int> b) {
	set<int> tmp;
	for(set<int>::iterator aitr = a.begin(); aitr != a.end(); aitr++) {
		if(b.find(*aitr) != b.end()) tmp.insert(*aitr);
	}
	return tmp.empty();
}

inline vector<ListDigraph::Node> intersectionTwoVectorsOfNodes(vector<ListDigraph::Node> a, vector<ListDigraph::Node> b) {
	vector<ListDigraph::Node> tmp(a.size());
	for(int i = 0; i < (int)a.size(); i++) {
		if(std::find(b.begin(), b.end(), a[i]) != b.end()) tmp.push_back(a[i]);
	}
	return tmp;
}

template <typename Element>
vector<Element> intersectionTwoVectors(vector<Element> a, vector<Element> b) {
	vector<Element> tmp; //(a.size());
	for(int i = 0; i < (int)a.size(); i++) {
		if(std::find(b.begin(), b.end(), a[i]) != b.end()) tmp.push_back(a[i]);
	}
	return tmp;
}


#ifdef __UOMAP__
inline set<int> getKeysFromMap(unordered_map<int, pair<int,int> > m) {
	set<int> keys;
	for(unordered_map<int, pair<int,int> >::iterator mitr = m.begin(); mitr != m.end(); mitr++) {
#else
inline set<int> getKeysFromMap(map<int, pair<int,int> > m) {
	set<int> keys;
	for(map<int, pair<int,int> >::iterator mitr = m.begin(); mitr != m.end(); mitr++) {
#endif
		keys.insert(mitr->first);
	}
	return keys;
}

//inline set<int> getKeysFromMap(unordered_map<int, pair<int,int> > m) {
//	set<int> keys;
//	for(unordered_map<int, pair<int,int> >::iterator mitr = m.begin(); mitr != m.end(); mitr++) {
//		keys.insert(mitr->first);
//	}
//	return keys;
//}

inline vector<int> getKeysFromMapAsVector(map<int, pair<int,int> > m) {
	vector<int> keys(m.size(), -1);
	int i = 0;
	for(map<int, pair<int,int> >::iterator mitr = m.begin(); mitr != m.end(); mitr++) {
		//keys.push_back(mitr->first);
		keys[i] = mitr->first;
		i++;
	}
	return keys;
}

inline int getIndexOfMaxValue(vector<int> v) {
	int maxVal = v[0];
	int maxInd = 0;
	for(int i = 0; i < (int)v.size(); i++) { if(v[i] > maxVal) { maxVal = v[i]; maxInd = i; } }
	return maxInd;
}

inline bool sortPair(pair<ListDigraph::Node, int> p1, pair<ListDigraph::Node, int> p2) { return p2.second > p1.second; }
inline bool sortPairInt(pair<int, int> p1, pair<int, int> p2) { return p2.first > p1.first; }
template <class part1, class part2>
inline bool sortPairGeneric(pair<part1, part2> p1, pair<part1, part2> p2) { return p2.first > p1.first; }
inline bool less_vectors(const vector<int>& a,const vector<int>& b) { return std::less<size_t>()(a.size(), b.size()); }
/**
 * is p1 a subpath of p2, ie-are all arcs of p1 contained in p2?
 */
inline bool isSubPath(vector<ListDigraph::Arc> p1, vector<ListDigraph::Arc> p2) {
	vector<ListDigraph::Arc> intArcs = intersectionTwoVectors(p1, p2);
	return (intArcs.size() == p1.size());
	//return intArcs == p1;
}

// for sake of time, taken from: http://codereview.stackexchange.com/questions/10130/edit-distance-between-two-strings
const size_t insert_cost = 1;
const size_t delete_cost = 1;
const size_t replace_cost = 1;

inline size_t min(size_t x, size_t y, size_t z)
{
    return x < y ? min(x,z) : min(y,z);
}

inline size_t edit_distance(const string& A, const string& B)
{
    size_t NA = A.size();
    size_t NB = B.size();

    vector<size_t> M0(NB+1), M1(NB+1);

    for (size_t b = 0; b <= NB; ++b)
        M0[b] = b * delete_cost;

    for (size_t a = 1; a <= NA; ++a)
    {
        M1[0] = a * insert_cost;

        for (size_t b = 1; b <= NB; ++b)
        {
            size_t x = M0[b] + insert_cost;
            size_t y = M1[b-1] + delete_cost;
            size_t z = M0[b-1] + (A[a-1] == B[b-1] ? 0 : replace_cost);
            M1[b] = min(x,y,z);
        }

        swap(M0,M1);
    }

    return M0[NB];
}

inline size_t delMut_distance(const string& A, const string& B)
{
    size_t NA = A.size();
    size_t NB = B.size();

    vector<size_t> M0(NB+1), M1(NB+1);

    for (size_t b = 0; b <= NB; ++b)
        M0[b] = b * delete_cost;

    for (size_t a = 1; a <= NA; ++a)
    {
        M1[0] = a * insert_cost;

        for (size_t b = 1; b <= NB; ++b)
        {
            size_t x = M0[b] + 0;	// remove any insert cost
            size_t y = M1[b-1] + delete_cost;
            size_t z = M0[b-1] + (A[a-1] == B[b-1] ? 0 : replace_cost);
            M1[b] = min(x,y,z);
        }

        swap(M0,M1);
    }

    return M0[NB];
}


#endif /* UTILS_H_ */
