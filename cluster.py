#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 19:08:52 2018

@author: dororo
"""

import numpy as np
from scipy.spatial import distance


class Edge(object):
    def __init__(self, b, d):
        self.b = b
        self.d = d

    def __str__(self):
        return "{}:{}".format(self.b, self.d)


class Point:
    def __init__(self, a):
        self.a = a
        self.adj = []

    def append(self, b, d):
        self.adj.append(Edge(b, d))

    def remove(self, index):
        rei = -1
        for i, e in enumerate(self.adj):
            if e.b == index:
                rei = i
            elif e.b >= index:
                e.b -= 1
        if rei >= 0:
            self.adj.pop(rei)

    def __getitem__(self, index):
        for e in self.adj:
            if e.b == index:
                return e

    def __str__(self):
        string = "{} -> ".format(self.a)
        for i, e in enumerate(self.adj):
            if i:
                string += " / "
            string += e.__str__()
        return string

    def __len__(self):
        return len(self.adj)


class Adjlist:
    def __init__(self):
        self.adjlist = []

    def append(self, item):
        self.adjlist.append(item)

    def pop(self, index):
        for adj in self.adjlist:
            if adj.a > index:
                adj.a -= 1
            adj.remove(index)
            if len(adj) == 0:
                self.adjlist.remove(adj)
        if index < len(self.adjlist):
            self.adjlist.pop(index)

    def __getitem__(self, index):
        return self.adjlist[index]

    def __str__(self):
        string = ""
        for i, adj in enumerate(self.adjlist):
            if i:
                string += "\n"
            string += adj.__str__()
        return string

    def __len__(self):
        return len(self.adjlist)+1


class Group:
    def __init__(self, n):
        self.g = [n]

    def __iadd__(self, other):
        self.g += other.g

    def append(self, other):
        for i in other.g:
            self.g.append(i)

    def __getitem__(self, index):
        return self.g[index]

    def __len__(self):
        return len(self.g)

    def __str__(self):
        string = ""
        for i in self.g:
            string += "{},".format(i)
        return string


class Cluster:
    def __init__(self, n):
        self.c = []
        for i in range(n):
            self.c.append(Group(i))

    def merge(self, a, b):
        self.c[a].append(self.c[b])
        self.c.pop(b)

    def __len__(self):
        return len(self.c)

    def __getitem__(self, index):
        return self.c[index]

    def display(self):
        for g in self.c:
            print(g)
        print("group size =", len(self))


def cossim(s1, s2, model, num_features, index2word_set):
    def avg_feature_vector(sentence):
        words = sentence.split()
        feature_vec = np.zeros((num_features, ), dtype='float32')
        n_words = 0
        for word in words:
            if word in index2word_set:
                n_words += 1
                feature_vec = np.add(feature_vec, model.wv[word])
        if n_words > 0:
            feature_vec = np.divide(feature_vec, n_words)
        return feature_vec
    s1_afv = avg_feature_vector(s1)
#    print("s1afv", s1_afv[0])
    s2_afv = avg_feature_vector(s2)
    return 1 - distance.cosine(s1_afv, s2_afv)


def eucsim(s1, s2, model, num_features, index2word_set):
    def avg_feature_vector(sentence):
        words = sentence.split()
        feature_vec = np.zeros((num_features, ), dtype='float32')
        n_words = 0
        for word in words:
            if word in index2word_set:
                n_words += 1
                feature_vec = np.add(feature_vec, model.wv[word])
        if n_words > 0:
            feature_vec = np.divide(feature_vec, n_words)
        return feature_vec
    s1_afv = avg_feature_vector(s1)
#    print("s1afv", s1_afv[0])
    s2_afv = avg_feature_vector(s2)
    return distance.euclidean(s1_afv, s2_afv)


def absdis(a, b):
    return abs(a-b)
