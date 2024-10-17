#include <omp.h>
#include <time.h>
#include <iostream>
#include <stack>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <set>
#include <vector>
#include <queue>
#include <map>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <bitset>
#include <math.h>

#define MAXINT ((unsigned) 4294967295)

using namespace std;

int dmax = 1, dmin = 1; // the number of two models， 初始测试要求和<=4
int totalV = 0, weig = 10, threads = 20;
int max1 = 20, max2 = 8, max3 = 7, max4 = 6;
int src, dst, aplace, bplace;
int dis = 1;
int maxdis = 0;
long long paths = 0;

double t1=0, t2=0, t3=0;


struct elem{

    u_int16_t a, b, c, d; // a: 下限值最大，剩余：和最小
    u_int16_t pdis;

    elem(){ a = 0, b = 0, c = 0, d = 0, pdis = 0; }

    elem(u_int16_t a1, u_int16_t b1, u_int16_t c1, u_int16_t d1){ // 固定写入
        a = a1, b = b1, c = c1, d = d1;
    }
    
    elem(u_int16_t a1, u_int16_t b1, u_int16_t c1, u_int16_t d1, u_int16_t dd){ // 固定写入
        a = a1, b = b1, c = c1, d = d1, pdis = dd;
    }

    void setvalue(u_int16_t val, int pla){ 
        switch (pla) {
            case 0: a=val; break;
            case 1: b=val; break;
            case 2: c=val; break;
            case 3: d=val; break;
            default: break;
        }
    }

    void setdis(u_int16_t val){
        pdis = val;
    }

    bool operator < (const elem& e1) const{
        u_int16_t val, val1;
        for (int i=0; i<dmax+dmin; ++i){
            switch (i) {
                case 0: val = a, val1 = e1.a; break;
                case 1: val = b, val1 = e1.b; break;
                case 2: val = c, val1 = e1.c; break;
                case 3: val = d, val1 = e1.d; break;
                default: break;
            }
            
            if (i<dmax){
                if (val < val1) return true;
                else if (val > val1) return false;
                else {}
            }else{
                if (val > val1) return true;
                else if (val < val1) return false;
                else {}
            }
        }

        return false;
	}

    bool operator <= (const elem& e1) const{
        // A <= B means that A is dominated by B
        u_int16_t val, val1;
        for (int i=0; i<dmax+dmin; ++i){
            switch (i) {
            case 0: val = a, val1 = e1.a; break;
            case 1: val = b, val1 = e1.b; break;
            case 2: val = c, val1 = e1.c; break;
            case 3: val = d, val1 = e1.d; break;
            default: break;
            }
            
            if (i<dmax){
                if (val > val1) return false;
            }else{
                if (val < val1) return false;
            }
        }

        return true;
	}

    bool operator <<= (const elem& e1) const{
        // a special case of A is dominated by B
        u_int16_t val, val1;
        for (int i=dmax; i<dmax+dmin; ++i){
            switch (i) {
            case 0: val = a, val1 = e1.a; break;
            case 1: val = b, val1 = e1.b; break;
            case 2: val = c, val1 = e1.c; break;
            case 3: val = d, val1 = e1.d; break;
            default: break;
            }
            
            if (val > val1) return false;
        }

        return true;
	}


    bool operator == (const elem& e1) const{
        u_int16_t val, val1;
        for (int i=0; i<dmax+dmin; ++i){
            switch (i) {
            case 0: val = a, val1 = e1.a; break;
            case 1: val = b, val1 = e1.b; break;
            case 2: val = c, val1 = e1.c; break;
            case 3: val = d, val1 = e1.d; break;
            default: break;
            }
            
            if (val != val1) return false;
        }

        return true;
	}
};


bool cmp(elem e1, elem e2){
    u_int16_t val1, val2;

    for (int i=0; i<dmax+dmin; ++i){
        switch (i) {
        case 0: val1 = e1.a, val2 = e2.a; break;
        case 1: val1 = e1.b, val2 = e2.b; break;
        case 2: val1 = e1.c, val2 = e2.c; break;
        case 3: val1 = e1.d, val2 = e2.d; break;
        default: break;
        }
        
        if (i<dmax){
            if (val1 > val2) return true;
            else if (val1 < val2) return false;
            else { }
        }else{
            if (val1 < val2) return true;
            else if (val1 > val2) return false;
            else { }
        }
    }

    return false;
}


bool cmpdis(elem e1, elem e2){

    if (e1.pdis > e2.pdis)
        return true;
    else if(e1.pdis < e2.pdis)
        return false;
    else{
        u_int16_t val1, val2;

        for (int i=0; i<dmax+dmin; ++i){
            switch (i) {
            case 0: val1 = e1.a, val2 = e2.a; break;
            case 1: val1 = e1.b, val2 = e2.b; break;
            case 2: val1 = e1.c, val2 = e2.c; break;
            case 3: val1 = e1.d, val2 = e2.d; break;
            default: break;
            }
            
            if (i<dmax){
                if (val1 > val2) return true;
                else if (val1 < val2) return false;
                else { }
            }else{
                if (val1 < val2) return true;
                else if (val1 > val2) return false;
                else { }
            }
        }
    }
    

    return false;
}


elem Merge(elem e1, elem e2){
    u_int16_t val1, val2;

    for (int i=0; i<dmax+dmin; ++i){
        switch (i) {
            case 0: {
                if (i<dmax) e1.a = min(e1.a, e2.a);
                else        e1.a += e2.a;
                break;
            }
            case 1:{
                if (i<dmax) e1.b = min(e1.b, e2.b);
                else        e1.b += e2.b;
                break;
            }
            case 2: {
                if (i<dmax) e1.c = min(e1.c, e2.c);
                else        e1.c += e2.c;
                break;
            }
            case 3: {
                if (i<dmax) e1.d = min(e1.d, e2.d);
                else        e1.d += e2.d;
                break;
            }
            default: break;
        }
    }
    

    return e1;
}


int hashelem(elem e1){
    return (int)(e1.a+e1.b*100+e1.c*100000+e1.d*100000000);
}


elem reverse_hash(int val){
    u_int16_t a, b, c, d;
    d = val / 100000000; 
    c = (val - d*100000000) / 100000;
    b = (val - d*100000000 - c*100000) / 100;
    a = val - d*100000000 - c*100000 - b*100;

    elem e(a, b, c, d);

    return e;
}


vector<elem> TSky; // 总的skyline值

vector<vector<elem> > label, rlabel;
vector<unordered_map<int, vector<elem> > > label_extra, rlabel_extra;
// vector<vector<pair<elem, elem> > > newLabels; // 两者维度均为 totalV，

vector<unordered_map<int, vector<int> > > newLabels; // 两者维度均为 totalV，

vector<map<int, vector<pair<unsigned, elem> > > > conR; // sketch graph

vector<vector<pair<unsigned, elem> > > con; // id + edge label
vector<pair<int, int> > v2degree;
vector<int> v2p, sList;

map<elem, int> Sky2Num;

void Print(elem e1){
    cout<<"("<<(int)(e1.a)<<", "<<(int)(e1.b)<<", "<<(int)(e1.c)<<", "<<(int)(e1.d)<<")  ";
    // cout<<"$["<<(int)(e1.b)<<", "<<(int)(e1.c)<<"]$"<<endl;
// <<(int)(e1.a)<<", "
}


void Print(vector<unsigned> vec){
    for (int i=0; i<vec.size(); ++i){
        cout<<vec[i]<<" ";
    }
    // cout<<endl;
}


void Print(vector<int> vec){
    for (int i=0; i<vec.size(); ++i){
        cout<<vec[i]<<" ";
    }
    // cout<<endl;
}


void ParaInitial(){
    label.resize(totalV);
    rlabel.resize(totalV);
    label_extra.resize(totalV);
    rlabel_extra.resize(totalV);
}


void SketchGraphPrint(){
    for (int i=0; i<totalV; ++i){
    
        map<int, std::vector<pair<unsigned, elem> >>& adj = conR[i];
        
        if (adj.size() == 0) continue;

        cout<<"id: "<<i<<endl;
        for (auto it=adj.begin(); it!=adj.end(); ++it){
            elem ee = reverse_hash(it->first);
            Print (ee);
            // Print (it->second);
        }
        cout<<endl;

    }
}


void GraphInitial(string filename){
    string s;
    const char *filepath = filename.c_str();
    ifstream infile;

    infile.open(filepath);
    if(!infile.is_open()){
        cout<<"No such file!"<<endl;
        exit(-1);
    }

    long xx = 0;
    srand((unsigned)time(NULL));
    
    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        char* s1 = strtok(strc," ");
        
        if (xx == 0){
            totalV = atoi(s1);
            con.resize(totalV);
            ParaInitial();
        }else{
            while(s1){
                unsigned va = xx - 1, vb = atoi(s1) - 1;
                u_int16_t a = dmax+dmin >= 1? u_int16_t((va+vb) % max1) + 1:0, 
                         b = dmax+dmin >= 2? u_int16_t((va+vb+1) % max2) + 1:0,
                         c = dmax+dmin >= 3? u_int16_t((va+vb+2) % max3) + 1:0, 
                         d = dmax+dmin >= 4? u_int16_t((va+vb+3) % max4) + 1:0;
                elem e(a, b, c, d); 
                con[va].push_back(make_pair(vb, e));

                // if (vb > va){
                //     cout<<" va: "<<va<<"   vb: "<<vb<<endl;
                //     Print(e);
                // }
                s1=strtok(NULL," ");
            }
        }

        xx += 1;
        
        delete s1, strc;
    }

    infile.close();
    long long edgs = 0;

    for( int i = 0; i < totalV; ++i ){
		v2degree.push_back(make_pair(con[i].size(), i)); // 
        edgs += con[i].size();
    }

	sort(v2degree.rbegin(), v2degree.rend());
    src = v2degree[aplace].second,  dst = v2degree[bplace].second;
    double cntt = (double)edgs / totalV;
    cout<<v2degree[0].first<<"  "<< cntt<<endl;
}


void Graphtest(string filename, string edgename){
    string s;
    const char *filepath = filename.c_str();
    ifstream infile;

    infile.open(filepath);
    if(!infile.is_open()){
        cout<<"No such file!"<<endl;
        exit(-1);
    }

    long xx = 0;
    srand((unsigned)time(NULL));
    
    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        char* s1 = strtok(strc," ");
        
        if (xx == 0){
            totalV = atoi(s1);
            con.resize(totalV);
            ParaInitial();
        }else{
            while(s1){
                unsigned va = xx - 1, vb = atoi(s1);
                // cout<<va<<"  "<<vb<<endl;
                elem e; 
                con[va].push_back(make_pair(vb, e));

                s1=strtok(NULL," ");
            }
        }

        xx += 1;
        
        delete s1, strc;
    }

    infile.close();

    const char *filepath1 = edgename.c_str();
    infile.open(filepath1);


    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        char* s1 = strtok(strc," ");
        unsigned va, vb, a=0, b=0, c=0, d = 0;
        xx = 0;

        while(s1){
            switch (xx) {
                case 0: va=atoi(s1); break;
                case 1: vb=atoi(s1); break;
                case 2: a=atoi(s1); break;
                case 3: b=atoi(s1); break;
                case 4: c=atoi(s1); break;
                case 5: d=atoi(s1); break;
                default: break;
            }

            xx += 1;
            s1=strtok(NULL," ");
        }

        va += 1;
        vb += 1;
        // cout<<va<<" "<<vb<<endl;
        for (int kk=0; kk<con[va-1].size(); ++kk){
            if (con[va-1][kk].first == vb-1){
                con[va-1][kk].second.a = a;
                con[va-1][kk].second.b = b;
                con[va-1][kk].second.c = c;
                con[va-1][kk].second.d = d;
            }
        }

        for (int kk=0; kk<con[vb-1].size(); ++kk){
            if (con[vb-1][kk].first == va-1){
                con[vb-1][kk].second.a = a;
                con[vb-1][kk].second.b = b;
                con[vb-1][kk].second.c = c;
                con[vb-1][kk].second.d = d;
            }
        }
        // cout<<a<<"  "<<b<<"  "<<c<<endl;
        
        delete s1, strc;
    }

    infile.close();
}


void CheckExist(vector<elem>& ins_Lab, vector<elem>& exist_Lab, 
           unordered_map<int, vector<elem> >& extra){
    int p = 0;

    for (int ii=0; ii<ins_Lab.size(); ++ii){
        elem& e_ins = ins_Lab[ii];
        int f1 = 0;
        elem edom = e_ins;
        
        for (int jj=0; jj<exist_Lab.size(); ++jj){
            elem& e_exist = exist_Lab[jj];

            if (e_ins == e_exist){ // extra应该也是放在一起的
                f1 = 1;
            }else if (e_ins <= e_exist){ 
                f1 = 1;
                if (e_ins <<= e_exist){
                    f1 = -1;
                    if (edom.a < e_exist.a)
                        edom.a = e_exist.a; // 找最优的
                }    
            }

            if (f1 == 1) break;
        }

        if (f1 == 0) 
            ins_Lab[p++] = ins_Lab[ii];

        if (f1 == -1){
            int hah = hashelem(edom), hah1 = hashelem(e_ins);
            extra[hah].push_back(e_ins);
            extra[hah].insert(extra[hah].end(), 
                              extra[hah1].begin(), extra[hah1].end());
        }
    }

    ins_Lab.resize(p);
}

void CheckSky(vector<elem>& ins_Lab, unordered_map<int, vector<elem> >& extra){
    int p = 0;

    for (int ii=0; ii<ins_Lab.size(); ++ii){
        elem& e_ins = ins_Lab[ii];
        int f1 = 0;
        
        for (int jj=0; jj<TSky.size(); ++jj){
            elem& e_exist = TSky[jj];

            if (e_ins <= e_exist or e_ins == e_exist){ // e_ins 是可以舍弃的
                f1 = 1; break;
            }
        }

        if (f1 == 0) 
            ins_Lab[p++] = ins_Lab[ii];
        else{
            int hah = hashelem(e_ins);
            extra[hah].clear();
            extra.erase(hah);
        }
    }

    ins_Lab.resize(p);
}


void RmvRedundant(vector<elem>& ins_Lab){
    if (ins_Lab.size() >= 2){
        int p = 1;
        for (int i=1; i<ins_Lab.size(); ++i){
            if (ins_Lab[i-1] == ins_Lab[i]) 
                continue;
            ins_Lab[p++] = ins_Lab[i];
        }
        ins_Lab.resize(p);
    }
}


void ResetSub(vector<elem>& Elements, unordered_map<int, vector<elem> >& extra){
    int p = 0, flg;
    for (int i=0; i<Elements.size(); ++i){
        elem& e1 = Elements[i];
        elem edom;
        flg = 0;

        for (int j=0; j<p; ++j){
            elem& e2 = Elements[j];

            if (e1 <= e2){
                flg = 1;
                // if (e1 <<= e2) {
                //     flg = -1, edom = e2;
                // } // 排序完成之后，就不用全部遍历了
                break; 
            }
        }

        if (flg == 0) Elements[p++] = Elements[i];
        if (flg == -1){
            
            int hah = hashelem(edom), hah1 = hashelem(e1);
            extra[hah].push_back(e1);

            if (extra.find(hah1) != extra.end())
            extra[hah].insert(extra[hah].end(), 
                              extra[hah1].begin(), extra[hah1].end());
        }
    }

    Elements.resize(p);
    vector<elem>(Elements).swap(Elements);
}




void ResetLabel(vector<elem>& Elements, unordered_map<int, vector<elem> >& extra){
    int p = 0, flg;
    for (int i=0; i<Elements.size(); ++i){
        elem& e1 = Elements[i];
        elem edom;
        flg = 0;

        for (int j=0; j<p; ++j){
            elem& e2 = Elements[j];

            if (e1 <= e2){
                flg = 1;
                if (e1 <<= e2) {
                    flg = -1, edom = e2;
                } // 排序完成之后，就不用全部遍历了
                break; 
            }
        }

        if (flg == 0) Elements[p++] = Elements[i];
        if (flg == -1){
            int hah = hashelem(edom), hah1 = hashelem(e1);
            extra[hah].push_back(e1);
            extra[hah].insert(extra[hah].end(), 
                              extra[hah1].begin(), extra[hah1].end());
        }
    }

    Elements.resize(p);

    for (auto it=extra.begin(); it!=extra.end(); ++it){
        Elements.insert(Elements.end(), it->second.begin(), it->second.end());
    }
    vector<elem>(Elements).swap(Elements);
}




void ResetSky(vector<elem>& Elements){
    int p = 0, flg;
    for (int i=0; i<Elements.size(); ++i){
        elem& e1 = Elements[i];
        
        flg = 0;

        for (int j=0; j<p; ++j){
            elem& e2 = Elements[j];

            if (e1 <= e2){
                flg = 1;
                // if (e1 <<= e2) flg = 0; 
            }

            if (flg == 1) break;
        }

        if (flg == 0) Elements[p++] = Elements[i];
    }

    Elements.resize(p);
    vector<elem>(Elements).swap(Elements);
}


void Inst(vector<elem>& E1, vector<elem>& E2){
    int a1=0, j;
    for (int i=0; i<E1.size(); ++i){
        for (j=a1; j<E2.size(); ++j){
            if (E2[j] < E1[i]){
                a1 = j+1;
                E2.insert(E2.begin()+j, E1[i]);
                break;
            }
        }

        if (j==E2.size()) E2.push_back(E1[i]);
    }

    vector<elem>().swap(E1);
}


void DstUpdate(int dis){
    
    vector<elem> Elems;

    for (int i=0; i<con[dst].size(); ++i){
        unsigned vid = con[dst][i].first;
        elem lb = con[dst][i].second;

        vector<elem>& labs = label[vid];

        if (labs.size() == 0) continue;

        for (int j=labs.size()-1; j>=0; --j){
            
            if (labs[j].pdis < dis-1) break;

            elem val = Merge(labs[j], lb);
            val.setdis(dis);
            Elems.push_back(val);
        }
    }

    if (Elems.size() > 0){
        sort(Elems.begin(), Elems.end(), cmp);
        Inst(Elems, TSky);
        ResetSky(TSky);  
    }
}




void ParallelBuild(){
    
    omp_set_num_threads(threads);
    vector<int> Cnt1(threads), Cnt2(threads);
    vector<vector<elem> > label_new(totalV), rlabel_new(totalV);

    label[src].push_back(elem(255, 0, 0, 0, 0));
    rlabel[dst].push_back(elem(255, 0, 0, 0, 0));
    double tt1 = 0, tt2 = 0, tt3 = 0;


    for( long long cnt1 = 1, cnt2 = 1; ; ++dis ){
        
        cnt1 = 0, cnt2 = 0;
        
        DstUpdate(dis); 

        #	pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();
            unordered_map<int, int> VVlab;
            for( int u = pid; u < totalV; u += np ){

                if (u == src or u == dst) continue;

                for (int i=0; i<con[u].size(); ++i){
                    unsigned vid = con[u][i].first;
                    elem lb = con[u][i].second;

                    vector<elem>& labs = label[vid];

                    for (int j=labs.size()-1; j>=0; --j){
                        
                        if (labs[j].pdis < dis-1) break;

                        elem val = Merge(labs[j], lb);  
                        val.setdis(dis);
                        label_new[u].emplace_back(val);
                    }

                    vector<elem>& rlabs = rlabel[vid];

                    for (int j=rlabs.size()-1; j>=0; --j){

                        if (rlabs[j].pdis < dis-1) break;

                        elem val = Merge(rlabs[j], lb);
                        val.setdis(dis);
                        rlabel_new[u].emplace_back(val);
                    }
                }

                if (label_new[u].size() > 0){
                    sort(label_new[u].begin(), label_new[u].end(), cmp);

                    RmvRedundant(label_new[u]);  
                    ResetSub(label_new[u], label_extra[u]); // 相同hop数，内部去冗余计算

                    CheckExist(label_new[u], label[u], label_extra[u]); // 基于已有子路径去重

                    CheckSky(label_new[u], label_extra[u]); // 基于已有skyline路径去重
                }

                if (rlabel_new[u].size() > 0){
                    sort(rlabel_new[u].begin(), rlabel_new[u].end(), cmp);
                    RmvRedundant(rlabel_new[u]); 
                    ResetSub(rlabel_new[u], rlabel_extra[u]); 

                    CheckExist(rlabel_new[u], rlabel[u], rlabel_extra[u]);
                    CheckSky(rlabel_new[u], rlabel_extra[u]); 

                }

                Cnt1[pid] +=  label_new[u].size();
                Cnt2[pid] += rlabel_new[u].size();
            }
        }

        #	pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();

            for ( int u = pid; u < totalV; u += np ){

                label[u].insert(label[u].end(), label_new[u].begin(), label_new[u].end());
                label_new[u].clear();
                vector<elem>().swap(label_new[u]);

                rlabel[u].insert(rlabel[u].end(), rlabel_new[u].begin(), rlabel_new[u].end());
                rlabel_new[u].clear();
                vector<elem>().swap(rlabel_new[u]);

            }
        }

        cnt1 = accumulate(Cnt1.begin(), Cnt1.end(), 0);
        cnt2 = accumulate(Cnt2.begin(), Cnt2.end(), 0);

        cout<<"dis: "<<dis<<"  "<<TSky.size()<<"  **  "<<cnt1<<"  "<<cnt2<<endl;
        
        Cnt1.clear(), Cnt1.resize(threads);
        Cnt2.clear(), Cnt2.resize(threads);



        if (cnt1 == 0 or cnt2 == 0) break;
    }


    #	pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for ( int u = pid; u < totalV; u += np ){
            sort(label[u].begin(), label[u].end(), cmp);
            ResetLabel(label[u], label_extra[u]);

            sort(rlabel[u].begin(), rlabel[u].end(), cmp);
            ResetLabel(rlabel[u], rlabel_extra[u]);
        }
    }
}



void LabelCombine(){
    
    newLabels.resize(totalV);

    sort(TSky.rbegin(), TSky.rend(), cmp);

    elem se(99, 0, 0, 0);
    int sval = hashelem(se);

    // cout<<endl<<"Skyline size: "<<TSky.size()<<endl;
    for (int i=0; i<TSky.size(); ++i){
        // Print(TSky[i]);
        // cout<<endl;
        
        Sky2Num[TSky[i]] = i+1;
        newLabels[src][sval].push_back(i+1);
        int dval = hashelem(TSky[i]);
        newLabels[dst][dval].push_back(i+1);

        if (TSky[i].pdis > maxdis)
            maxdis = TSky[i].pdis;
    }

    

    omp_set_num_threads(threads);
    #	pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for ( int u = pid; u < totalV; u += np ){

            if (u == src or u == dst) continue;
            
            vector<elem>& sLab = label[u];
            vector<elem>& dLab = rlabel[u];
 
            if (sLab.size() == 0 or dLab.size() == 0) continue;

            sort(sLab.rbegin(), sLab.rend(), cmp);
            sort(dLab.rbegin(), dLab.rend(), cmp);

            for (int ii=0; ii<sLab.size(); ++ii){
                
                elem& s1 = sLab[ii];

                for (int jj=0; jj<dLab.size(); ++jj){
                    
                    elem& d1 = dLab[jj];

                    elem newp = Merge(s1, d1);

                    for(int kk=0; kk<TSky.size(); ++kk){ // 
                        
                        elem ele = TSky[kk];

                        if (dmax > 0 and newp.a < ele.a) break;

                        if (dmax == 0 and newp.a > ele.a) break;

                        if (ele == newp){
                            int hashid = hashelem(s1);
                            newLabels[u][hashid].push_back(Sky2Num[ele]);
                            v2p[u] += 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    
}


void GraphReduction(){
    
    conR.resize(totalV);

    omp_set_num_threads(threads);

    # pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for( int u = pid; u < totalV; u += np ){

            if (newLabels[u].size() == 0) continue;

            vector<pair<unsigned, elem> >& adj = con[u];

            for (int i=0; i<adj.size(); ++i){

                unsigned v = adj[i].first;
                elem elab = adj[i].second;

                if (newLabels[v].size() == 0) continue;

                if (u == v) continue;

                for (auto it=newLabels[u].begin(); it!=newLabels[u].end(); ++it){
                    int hashval = it->first;
                    elem aa = reverse_hash(hashval);
                    elem newlab = Merge(elab, aa);
                    int hashVV = hashelem(newlab);

                    if (newLabels[v].find(hashVV) != newLabels[v].end()){
                        v2p[u] += 1;
                        conR[u][hashval].push_back(make_pair(v, elab)); // 满足访问条件
                        // cout<<u<<"  -  "<<v<<endl;
                    }
                }
            }
    
        }
    }

}



void DFSearch(unsigned u, unsigned tgt, 
              vector<unsigned>& s_, vector<elem>& lab_){
    
    s_.push_back(u);
    elem tag = lab_[lab_.size()-1];
    int hashid = hashelem(tag);

    vector<pair<unsigned, elem> >& adj = conR[u][hashid];

    if (s_.size() == maxdis - 1){
        paths += adj.size();
    }else{
        for (int i=0; i<adj.size(); ++i){

            unsigned v = adj[i].first;
            elem elab = adj[i].second;

            if (v == dst){
                paths += 1;
            }else{
                elem new_sub = Merge(tag, elab);
                lab_.push_back(new_sub);
                DFSearch(v, tgt, s_, lab_);
                lab_.pop_back();
            }
        }
    }

    s_.pop_back();
    
    return;
}







void test(){
    elem e1(3, 462, 0, 0);

    int aa = hashelem(e1);
    elem cc = reverse_hash(aa);
    Print(cc);
}




int main(){

    string na = "sk2005";
    string name = "/mnt/data/zyy/Full/"+na+".graph";
    string edgename = "/mnt/data/zyy/Full/"+na+"_edge.attribute";

    max1 = 10, max2 = 10, max3 = 10, max4 = 10;
    dmax = 1, dmin = 1;
    aplace = 1, bplace = 10;

    vector<unsigned> stk;
    vector<elem> lab;  
    lab.push_back(elem(99, 0, 0, 0));

    GraphInitial(name);
    
    // Graphtest(name, edgename);
    // src = 58718, dst = 766226;
    // // src = 92973,  dst = 4522357;
  
    // v2p.resize(totalV, 0);
    // cout<<src<<"  "<<dst<<endl;

    // threads = 40;

    // srand((unsigned)time(NULL));
    
    // double tt1 = omp_get_wtime();

    // ParallelBuild();
    // cout<<"Compute end!"<<endl;
    // double tt2 = omp_get_wtime();

    // LabelCombine();
    // // cout<<"Combine end!  "<<maxdis<<endl;
    // double tt3 = omp_get_wtime();
    // GraphReduction();
    
    // vector<vector<elem> >().swap(label), vector<vector<elem> >().swap(rlabel);
    // vector<unordered_map<int, vector<int> > >().swap(newLabels);

    // cout<<"Condense end!"<<endl;
    // double tt4 = omp_get_wtime();
    // int aa = accumulate(v2p.begin(), v2p.end(), 0);
    // cout<<"aa: "<<aa<<"   "<<maxdis<<endl;
    // DFSearch(src, dst, stk, lab);
    
    // cout<<tt2-tt1<<" "<<tt3-tt2<<"  "<<tt4-tt3<<endl;

    // cout<<"time: "<<omp_get_wtime()-tt1<<" s  ||  Paths:"<<paths<<endl;

    return 0;
}
