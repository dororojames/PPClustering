#include <stddef.h>
#include <stdlib.h>
#include "cluster.h"

struct edge {

    // pointers of point 
    Point* b;

    // the distance of edges
    double d;
};

struct point {
    // point id
    unsigned int id;
    
    // array of pointer of edges
    Edge** adj;

    // how many edges the point has
    unsigned int length;
};

struct group
{
    // group id
    unsigned int id;

    // array of pointer of points
    Point** p;

    // how many points in the group
    unsigned int length;
};

struct cluster {
    // array of pointer of groups
    Group** g;

    // how many group in the cluster 
    unsigned int length;

    // array of pointer of points, for easy access to the point
    Point** p;

    // how many points in the cluster 
    unsigned int amount_point;

    /*
        adjacent list, which store the distance
        ex. adj[0][1] is the distance between point 0 and point 1
    */
    double** adj;
};

/* 
    funct: create point and allocate the memory it needs
    parms:
        id: the id of point
    
    return:
        a pointer of Point
*/
Point* Create_Point(unsigned int id, unsigned int n){
    Point* p = (Point*) malloc(sizeof(Point*));
    p->id = id;
    p->adj = (Edge**) calloc(n, sizeof(Edge*));
    p->length = 0;

    return p;
}

/* 
    funct: create group and allocate the memory it needs
    parms:
        id: the id of group
    
    return:
        a pointer of Group
*/
Group* Create_Group(unsigned int id, unsigned int n){
    Group* g = (Group*) malloc(sizeof(Group*));
    g->id = id;
    g->p = (Point**) malloc(1 * sizeof(Point*));
    g->length = 1;
    
    g->p[0] = Create_Point(id, n);

    return g;
}


/* 
    funct: create cluster and allocate the memory it needs
    parms:
        n: the amount of points
    
    return:
        a pointer of Cluster
*/
Cluster* Create_Cluster(unsigned int n){
    Cluster* c = (Cluster*) malloc(sizeof(Cluster*));
    c->g = (Group**) malloc(n * sizeof(Group*));
    c->length = n;

    c->p = (Point**) malloc(n * sizeof(Point*));
    c->amount_point = n;
    c->adj = (double**) malloc(n * sizeof(double*));

    int i;
    for(i = 0; i < n; ++i){
        c->g[i] = Create_Group(i, n);
        c->p[i] = c->g[i]->p[0];

        c->adj[i] = malloc(n * sizeof(double));
    }

    return c;
}

/* 
    funct: append the edge to target point
    parms:
        p: pointer of target point
        b: pointer of point on the other end of the edge
        d: the distance of the edge
    return:
        none
*/
void Append(Point* p, Point* b, double d){
    Edge* e = (Edge*) malloc(sizeof(Edge*));
    e->b = b;
    e->d = d;

    p->adj[b->id] = e;
    ++p->length;
}

/*
    funct: insert the edge into the cluster, this function is supposed to be the initialization function,
           which means it is called when inserting the data to the cluster after the cluster is created. 
    parms:
        p1: pointer of starting point of the edge
        p2: pointer of end point of the edge
        d: the distance of the edge
    return:
        none
*/
void Inster_Edge(Cluster* c, unsigned int p1, unsigned int p2, double d){
    Append(c->p[p1], c->p[p2], d);

    c->adj[p1][p2] = d; 
}

/* 
    funct: merge two groups
    parms:
        self: the group that absorb the other group
        other: the group that is going to be absorbed
    return:
        none
*/
void Merge_Group(Group* self, Group* other){
    unsigned int insert_position = self->length;
    self->length += other->length;
    self->p = realloc(self->p, self->length * sizeof(Point*));

    int i;
    for(i = 0; i < other->length; ++i){
        self->p[insert_position + i] = other->p[i];
    }

    free(other->p);
    free(other);
}

/* 
    funct: merge two group in the cluster given by two group id
    parms:
        c: the cluster
        g1: group id of the group that will absorb the other group
        g2: group id of the group that will be absorbed
    return:
        none
*/
void Merge(Cluster* c, unsigned int g1, unsigned int g2){
    Merge_Group(c->g[g1], c->g[g2]);
    
    c->g[g2] = NULL;
}


void Clear(Cluster* c){
    int i, j, k;
    for(i = 0; i < c->length; ++i){
        if(c->g[i] == NULL){
            continue;
        }
        else{
            for(j = 0; j < c->g[i]->length; ++j){
                for(k = 0; k < c->g[i]->p[j]->length; ++k){
                    free(c->g[i]->p[j]->adj[k]);
                }
                
                free(c->g[i]->p[j]->adj);
                free(c->g[i]->p[j]);
            }

            free(c->g[i]);
        }
    }
    
    for(i = 0; i < c->amount_point; ++i){
        free(c->adj[i]);
    }
    free(c->adj);
}