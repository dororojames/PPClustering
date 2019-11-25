#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# encoding=utf-8
"""
Created on Sun May 13 14:53:50 2018

@author: dororo
"""

import copy
import os
import time

import jieba
import jieba.analyse
import jieba.posseg as jbps
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from gensim.models import Word2Vec
from opencc import OpenCC
from sklearn.decomposition import PCA
from textrank4zh import TextRank4Sentence

import cluster
import savetext
import selectionalgo

_PATH = os.path.dirname(os.path.abspath(__file__))+"/"
htmldir, savedir, tempdir, sentences, indexsentences = "", "", "", [], []


class Sentence(object):
    def __init__(self, sen, line, d, p, l, score):
        self.sentence = sen
        self.sentencelist = " ".join(line)
        self.doc = d
        self.para = p
        self.lineindex = l
        self.score = score

    def __str__(self):
        return "{} {} {} {} {}".format(self.sentence[:10], self.score, self.doc, self.para, self.lineindex)

    def __repr__(self):
        return self.__str__()


def dircheck(path, tstamp):
    if not os.path.isdir(path):
        os.mkdir(path)
    if tstamp != "":
        path += tstamp+"/"
        os.mkdir(path)
    return path


def getwebdata(search_name):
    t0 = time.time()
    sitelist = selectionalgo.selectalgo(search_name, _PATH)
    if sitelist == "noarticletext":
        return False
    print("search time", time.time()-t0)
    for file in os.listdir(savedir):
        os.rename(savedir+file, tempdir+file)
    for file in os.listdir(htmldir):
        os.rename(htmldir+file, tempdir+file)
    for site in sitelist:
        site.display()
        with open(htmldir+site.webname+".html", "w") as fp:
            fp.writelines(site.soup.prettify())
        sitetext = savetext.savetext(site.soup)
        with open(savedir+"{}_{}.txt".format(site.webname, site.score), "w") as outfp:
            outfp.writelines(sitetext)
    return True


def is_alphabet(uchar):
    """判断一个unicode是否是英文字母"""
    return (uchar >= u'\u0041' and uchar <= u'\u005a') or (uchar >= u'\u0061' and uchar <= u'\u007a')


def loadtext(textdir, mixedversion=True, textrankrate=0.5):
    t0 = time.time()
    docs = []
    filelist = os.listdir(textdir)
    for f in filelist:
        if not is_alphabet(f[0]):
            filelist.remove(f)
    # print(filelist)
    tr4s = TextRank4Sentence()
    for file in sorted(filelist, key=lambda file: float(file.split("_")[1][:-4]), reverse=True):
        print(file)
        with open(savedir+file) as txtfile:
            eofp = open(tempdir+file, "w")
            contents = txtfile.read()
            if mixedversion:
                tr4s.analyze(text=contents, lower=True, source='all_filters')
                docsum = tr4s.get_key_sentences(
                    num=int(len(tr4s.sentences)*textrankrate))
                # print(docsum)
                contents = [item.sentence for item in docsum]
            else:
                contents = contents.splitlines()
            # print(contents)
            paragraph, sentence = [], []
            for content in contents:
                line, i, l = "", 0, len(content)
                while i < l:
                    if content[i] == " " and i+1 < l and not is_alphabet(content[i+1]):
                        i += 1
                        continue
                    if content[i] == "　":
                        i += 1
                        continue
                    line += content[i]
                    if content[i] in ["。", "？", "！", "?", "!"]:
                        if i+1 < l and content[i+1] in ["」", "』", "”"]:
                            line += content[i+1]
                            i += 1
                        if line[0] in ["(", "（"] and sentence:
                            sentence[-1] += line
                            line = ""
                        elif i+1 < l and content[i+1] not in ["》", "〉", ")", "）"]:
                            sentence.append(line)
                            line = ""
                    i += 1
                if line:
                    sentence.append(line)
                if sentence:
                    for s in sentence:
                        eofp.write(s+"\n")
                    paragraph.append(sentence)
                    sentence = []
            eofp.close()
            docs.append(paragraph)
    print("txt loding completed", time.time()-t0)
    xlist = ["。", "？", "！", "～", "?", "!", " ", "　", "」", "』", "”",
             "(", "（", "》", "〉", ")", "）", "，", "：", "」", "、", "《",
             "；", "「", "%"]
    t0 = time.time()
    dic = {}
    with open(_PATH+"text segmentation/1998.csv", "r", encoding="utf-8") as dictxt:
        lines = dictxt.read().splitlines()
        for line in lines:
            line = line.split(",")
            dic[line[1]] = float(line[3])
    sentences = []
    indexsentences = []
    opcc = OpenCC('s2twp')
    for d, doc in enumerate(docs[:]):
        for p, para in enumerate(doc):
            for sid, sen in enumerate(para):
                # print(sen)
                sentence = []
                idf = count = w = score = 0
                hitset = set()
                words = jbps.cut(sen)
                for word, flag in words:
                    w += 1
                    word = opcc.convert(word)
                    try:
                        float(word)
                    except:
                        if not word in xlist:
                            # print(word, flag, end="|")
                            sentence.append(word)
                        if not word in dic or word in hitset:
                            continue
                        if flag in ["l", "n", "nr"] and dic[word] >= 30:
                            idf += 1/dic[word]
                            count += 1
#                print("\n", sentence, "\n")
                if w >= 10 and count > 3:
                    score = idf/count
                if sentence and score > 0:
                    sentences.append(sentence)
                    indexsentences.append(
                        Sentence(sen, sentence, d, p, sid, score))
    print("segmentation completed", time.time()-t0)
    print("total sentences:", len(sentences))
#    print(sentences)
    # for s in indexsentences:
    #     print(s)
    return sentences, indexsentences


def drawW2Vtrainningresult(model, w2vtotaltraintime, w2vdir, timestamp):
    words = list(model.wv.vocab)
    vector2D = PCA(n_components=2).fit_transform(model.wv.vectors)
    zhfont = matplotlib.font_manager.FontProperties(
        fname='/System/Library/Fonts/PingFang.ttc')
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111)
    xmax, xmin, ymax, ymin = 0, 0, 0, 0
#    ylimit = [-1.5, 1]
#    xlimit = [-20, 15]
    for i, v in enumerate(vector2D):
        word = words[i]
#        if not ylimit[0]<=v[1]<=ylimit[1] or not xlimit[0]<=v[0]<=xlimit[1]:
#            continue
        xmax = max(xmax, v[0])
        xmin = min(xmin, v[0])
        ymax = max(ymax, v[1])
        ymin = min(ymin, v[1])
        ax.text(v[0], v[1], word, fontproperties=zhfont)
    print(xmax, xmin, ymax, ymin)
    ax.axis([xmin, xmax*1.05, ymin, ymax*1.05])
    plt.savefig(w2vdir+"W2V2D_{}_{}.png".format(w2vtotaltraintime, timestamp))
#    plt.show()


def clustering(search_name, summarytype, trainning=False, timestamp=""):
    #    ----------word2vec train---------
    w2vdir = _PATH+"word2vec/{}/".format(search_name)
    if not os.path.isdir(w2vdir):
        os.mkdir(w2vdir)
        trainning = True
    w2vtotaltraintime, w2vdimsize = 200000, 100
    t0 = time.time()
    if trainning:
        print("begin word2vector trainning")
        model = Word2Vec(sentences, size=w2vdimsize, min_count=1,
                         workers=2, iter=w2vtotaltraintime)
        model.save(w2vdir+"word2vec_{}.model".format(w2vtotaltraintime))
        print("w2v trainning {} completed".format(
            w2vtotaltraintime), time.time()-t0)
        # -----------draw train result-----
        drawW2Vtrainningresult(model, w2vtotaltraintime, w2vdir, timestamp)
    else:
        # print("load model with {} train".format(w2vtotaltraintime))
        # model = Word2Vec.load(
        #     w2vdir+"word2vec_{}.model".format(w2vtotaltraintime))
        model = Word2Vec(sentences, size=w2vdimsize, min_count=1, workers=2)
        # -----------draw train result-----
        drawW2Vtrainningresult(model, w2vtotaltraintime, w2vdir, timestamp)
    print("totalwords=", len(model.wv.vectors))
#    ------------clustering-----------
    t0 = time.time()
    index2word_set = set(model.wv.index2word)
    adjlist = cluster.Adjlist()
    for i in range(len(indexsentences)-1):
        adjlist.append(cluster.Point(i))
        for j in range(i+1, len(indexsentences)):
            # print(indexsentences[i].sentence, indexsentences[j].sentence)
            d = cluster.cossim(indexsentences[i].sentencelist,
                               indexsentences[j].sentencelist,
                               model, w2vdimsize, index2word_set)
            # print(d)
            adjlist[i].append(j, d)
    print("distances calculated completed", time.time()-t0)
    print("begin clustering")
    csdir = dircheck(tempdir+"sc/", tstamp="")
    groupdir = dircheck(tempdir+"group/", tstamp="")
    # print(adjlist)
    clustersize = int(len(sentences)*0.45)
    opcc = OpenCC('s2twp')
    t00 = time.time()

    groupadj = copy.deepcopy(adjlist)
    clusters = cluster.Cluster(len(groupadj))
    while len(clusters) >= clustersize:
        # clusters.display()
        gl = len(clusters)
        t0 = time.time()
        mind, ga, gb = 9999, 0, 0
        for i in range(len(groupadj)-1):
            for j in range(i+1, len(groupadj)):
                gd = groupadj[i][j].d
                if gd < mind:
                    mind, ga, gb = gd, i, j
        # print(gl, "findtime", time.time()-t0, end=" ")
        # t0 = time.time()
        # print("merge{:4d}{:4d} mind= {}".format(ga, gb, mind))
        for i in range(len(clusters)):
            if i < ga:
                groupadj[i][ga].d = max(groupadj[i][ga].d, groupadj[i][gb].d)
            elif i > ga:
                for b in range(ga+1, gl):
                    if b == gb:
                        continue
                    ra, rb = min(gb, b), max(gb, b)
                    groupadj[ga][b].d = max(
                        groupadj[ga][b].d, groupadj[ra][rb].d)
                break
    #    print("change weighttime", time.time()-t0, end=" ")
    #    t0 = time.time()
        if gl == clustersize or (gl < len(sentences)//2 and (gl-clustersize) % 10 == 0):
            print(gl)
            cslist = []
            for i, g in enumerate(clusters):
                with open(groupdir+"{}g{:02d}.txt".format(gl, i), "w") as gfp:
                    for sid in g:
                        gfp.write(indexsentences[sid].sentence+"\n")
                    md, csid = 99999, 0
                    if len(g) == 2:
                        if indexsentences[g[0]].score < indexsentences[g[1]].score:
                            csid = g[1]
                        else:
                            csid = g[0]
                    else:
                        for lsid in g:
                            cd = 0
                            for rsid in g:
                                if lsid == rsid:
                                    continue
                                ia, ib = min(lsid, rsid), max(lsid, rsid)
                                cd += adjlist[ia][ib].d
                            if cd < md:
                                md, csid = cd, lsid
                    gfp.write(
                        "cs="+indexsentences[csid].sentence+" {}\n".format(indexsentences[csid].score))
                    cslist.append(indexsentences[csid])
            cslist = sorted(cslist, key=lambda s: s.score, reverse=True)[:20]
            for cs in cslist:
                if cs.score < 0.0075:
                    print(cs)
                    cslist.remove(cs)
            with open(csdir+"{}sc.txt".format(gl), "w") as csfp:
                for s in sorted(cslist, key=lambda s: (s.doc, s.para, s.lineindex)):
                    csfp.write(opcc.convert(s.sentence)+"\n")
            if gl == clustersize:
                with open(_PATH+"summary/"+search_name+"_summary"+summarytype+timestamp+".txt", "w") as csfp:
                    for s in sorted(cslist, key=lambda s: (s.doc, s.para, s.lineindex)):
                        csfp.write(opcc.convert(s.sentence)+"\n")
                break
        clusters.merge(ga, gb)
        groupadj.pop(gb)
        # print(gl, "mergetime", time.time()-t0)
        # print(groupadj)
    # clusters.display()
    print("clustering completed", time.time()-t00)


def getsummary(search_name="人工智慧", isgetwebdata=True, mixver=False, trrate=0.7, istrainning=True):
    global savedir, htmldir, tempdir, sentences, indexsentences
    tstamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    savedir = dircheck(_PATH+"savetext/{}/".format(search_name), tstamp="")
    htmldir = dircheck(_PATH+"html/{}/".format(search_name), tstamp="")
    tempdir = dircheck(_PATH+"temp/{}/".format(search_name), tstamp="")
    print("開始搜尋："+search_name)
    stype = ""
    if mixver:
        stype += "_mix_{}_".format(trrate)
    else:
        stype += "_cluster_"
    print("summary type is", stype)
    if istrainning:
        print("trainning")
    else:
        print("no trainning")
    # ----------getwebdata-----------
    if isgetwebdata:
        wikiresult = getwebdata(search_name)
        if not wikiresult:
            return None
    # ----cut sentence and word------
    sentences, indexsentences = loadtext(
        savedir, mixedversion=mixver, textrankrate=trrate)
    # -------------clustering--------
    clustering(search_name=search_name, summarytype=stype,
               trainning=istrainning, timestamp=tstamp)


if __name__ == "__main__":
    search_name = "人工智慧"
    #    search_name = input("輸入搜尋項目名稱:")

#    jieba.set_dictionary(_PATH+"text segmentation/dict.txt.big")
    mixver, trrate, istrainning = False, 0.7, True
    mixverinput = input("ismix?(y/N)")
    if str.lower(mixverinput) == "y":
        mixver = True
        trrateinput = input("Please input trrate(0.3~0.7)>>>")
        trrate = min(max(float(trrateinput), 0.3), 0.7)
    traininput = input("istrainning?(Y/n)")
    if str.lower(traininput) == "n":
        istrainning = False
    getsummary(search_name=search_name, isgetwebdata=True,
               mixver=mixver, trrate=trrate, istrainning=istrainning)
