#ifndef __Cluster__
#define __Cluster__
typedef struct edge Edge;
typedef struct point Point;
typedef struct group Group;
typedef struct cluster Cluster;

extern Point* Create_Point(unsigned int id, unsigned int n);
extern Group* Create_Group(unsigned int id, unsigned int n);
extern Cluster* Create_Cluster(unsigned int n);

extern void Append(Point* p, Point* b, double d);
extern void Inster_Edge(Cluster* c, unsigned int p1, unsigned int p2, double d);
extern void Merge_Group(Group* self, Group* other);
extern void Merge(Cluster* c, unsigned int g1, unsigned int g2);
extern void Clear(Cluster* c);
#endif // __Cluster__