﻿#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <mutex>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
#include <algorithm>
using namespace std;

unordered_map<string, string> Mimazi;//密码子表
string referenceSeq = "+";//参考基因组
string referenceRBD_Nuc = "";//参考rbd 331 531 氨基酸组成
string referenceRBD_AA = "";//参考rbd 331 531 氨基酸组成

//string split函数
vector<string> splitStr(string str, char delimiter)
{
    vector<string> r;
    string tmpstr;
    int i, j, k;
    vector<int> pointList;
    pointList.emplace_back(-1);
    for (i = 0; i < str.length(); i++)
        if (str[i] == delimiter)pointList.emplace_back(i);
    pointList.emplace_back(i);
    for (i = 1; i < pointList.size(); i++)
        r.emplace_back(str.substr(pointList[i - 1] + 1, (pointList[i] - pointList[i - 1] - 1)));
    return r;
}

void ReadinMimazi(string fold)//readin codon and reference genome 读取密码子和参考基因组
{
    ifstream read(fold + "/mimazi.txt");
    if (!read.is_open())
    {
        printf("/// ERROR Mut2Haplo: Can NOT open mimazi.txt ///");
        return;
    }
    string line;
    while (getline(read,line))
    {
        vector<string> line1 = splitStr(line, '\t');
        Mimazi.insert(pair<string,string>(line1[0], line1[1]));
    }
    read.close();

    //参考基因组 reference genome
    ifstream read2(fold + "/reference.fa");
    getline(read2, line);
    getline(read2, line);
    referenceSeq += line;
    read2.close();

    ifstream read3(fold + "/refRBD.fa");
    getline(read3, line);
    getline(read3, line);
    referenceRBD_AA += line;
    read3.close();

    referenceRBD_Nuc = referenceSeq.substr(22553, 23155 - 22553 + 1);
    return;
}

string Mutation2AA_RBD(string Mut)//输入核酸突变，返回氨基酸突变,限定在RBD，其它位置的突变返回“NA” transfer nuc mut 2 AA mut
{
    int i, j, k;
    vector<string> line2 = splitStr(Mut, '/');
    k = line2[0].length() - 1;
    while (!(line2[0][k] <= '9') || !(line2[0][k] >= '0')) k--;
    string refNuc = line2[0].substr(k + 1, line2[0].length() - k - 1);
    string mutNuc = line2[1];
    int pos = stoi(line2[0].substr(0, k + 1));

    if (refNuc.length() + mutNuc.length() != 2)//ignore inDels
        return "NA";

    string ATCG = "ATCG";

    if (pos >= 22553 && pos <= 23155 && ATCG.find(mutNuc)!=-1)//choose positions
    {
        int forward = (pos - 22517) % 3;
        int backward = 2 - forward;
        string codonRef = "";
        string codonMut = "";
        for (i = pos - forward; i < pos; i++)
        {
            codonRef += referenceSeq[i];
            codonMut += referenceSeq[i];
        }
        codonRef += refNuc;
        codonMut += mutNuc;
        for (i = pos + 1; i <= pos + backward; i++)
        {
            codonRef += referenceSeq[i];
            codonMut += referenceSeq[i];
        }

        string refAmino = Mimazi[codonRef];
        string mutAmino = Mimazi[codonMut];
        if (refAmino != mutAmino)
        {
            string aminoChange = refAmino + to_string((pos - 21563) / 3 + 1) + mutAmino;
            return aminoChange;
        }
        else
            return "NA";
    }
    else
        return "NA";
}

string MutationList2AA_RBD(string Mutlist)////输入核酸突变，返回氨基酸突变,限定在RBD，纳入考虑单个密码子中发生多个位点突变的情况
{
    int i, j, k;
    string transfered = "";
    string nucSeq = referenceRBD_Nuc;
    vector<string> Mut = splitStr(Mutlist, ',');
    for (j = 0; j < Mut.size(); j++)
    {
        vector<string> line2 = splitStr(Mut[j], '/');
        k = line2[0].length() - 1;
        while (!(line2[0][k] <= '9') || !(line2[0][k] >= '0')) k--;
        string refNuc = line2[0].substr(k + 1, line2[0].length() - k - 1);
        string mutNuc = line2[1];
        int pos = stoi(line2[0].substr(0, k + 1));

        if (refNuc.length() + mutNuc.length() != 2)//ignore inDels
            continue;

        string ATCG = "ATCG";

        if (pos >= 22553 && pos <= 23155 && ATCG.find(mutNuc) != -1)//choose positions
        {
            nucSeq[pos - 22553] = mutNuc[0];
        }
    }
    for (i = 0; i < nucSeq.size(); i += 3)
    {
        string codonRef = "" + referenceRBD_Nuc.substr(i,1) + referenceRBD_Nuc.substr(i+1, 1) + referenceRBD_Nuc.substr(i+2, 1);
        string codonMut = "" + nucSeq.substr(i, 1) + nucSeq.substr(i + 1, 1) + nucSeq.substr(i + 2, 1);
        if (Mimazi[codonRef] != Mimazi[codonMut])
        {
            transfered += Mimazi[codonRef] + to_string(i / 3 + 331) + Mimazi[codonMut] + " ";
        }
    }
    if (transfered == "")
        return "";
    else
        return transfered.substr(0, transfered.size() - 1);
}

void ConvertAndOutput(string file)
{
    int i, j, k;
    string line;
    ifstream read(file);
    if (!read.is_open())
    {
        printf("/// ERROR Mut2Haplo: Can NOT open the MutMetadata file ///");
        return;
    }
    ofstream write(file + ".RBDAAmut");
    if (!write.is_open())
    {
        printf("/// ERROR Mut2Haplo: Can NOT create the Output file ///");
        return;
    }
    getline(read, line);
    while (getline(read,line))
    {
        vector<string> line1 = splitStr(line, '\t');
        vector<string> mut1 = splitStr(line1[1], ' ');
        string mutation = "";
        for (i = 0; i < mut1.size(); i++)
        {
            string covmut = Mutation2AA_RBD(mut1[i]);
            if (covmut != "NA")
                mutation += covmut + " ";
        }
        if (mutation.size() > 1)
            mutation = mutation.substr(0, mutation.size() - 1);
        string output = line1[0] + "\t";
        output += mutation + "\t";
        output += line1[2] + "\t";
        output += line1[3] + "\t";
        output += line1[4];
        write << output << endl;
    }
    return;
}

int main(int argc, char* argv[])//1个输入，需要转换的目标文件; 
{
    if (argc != 2)
    {
        printf("/// ERROR Mut2Haplo: 1 input needed. The path of MutMetadata file ///");
        return 1;
    }
    string fold = argv[0];
    int k = fold.size();
    while (fold[k] != '/') k--;
    fold = fold.substr(0, k);

    ReadinMimazi(fold);
    ConvertAndOutput(argv[1]);
    //ReadinMimazi("/mnt/g/VariationMutation/VirusAntibodyEscape/LinuxPipeline/SelectiveCoEfficient/Script");
    //ConvertAndOutput("/mnt/g/VariationMutation/回复突变/iqtreePopulation20221029/Linux/sample_mut_loc_time.tsv");

    return 0;
}