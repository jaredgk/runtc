#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>


using namespace std;



vector<string> split(string s, string delim) {
    vector<string> out;
    if(s.size() == 0) { return out; }
    int pos = 0;
    string token;
    while ((pos = s.find(delim)) != string::npos) {
        token = s.substr(0,pos);
        out.push_back(token);
        s.erase(0,pos+delim.length());
    }
    out.push_back(s);
    return out;
}

class snp_data {
    public:
        int allele_count;
        int position;
        vector<char> alleles;
        snp_data(int ac, int pos, vector<char> as) {
            allele_count = ac;
            position = pos;
            alleles = as;
        }
};

double getGenCol(string line, int col, double & phys) {
    stringstream s(line);
    s >> phys;
    double j = 0, out = 0;
    for(int i = 0; i < col; i++) {
        s >> j;
    }
    bool b = (s >> out);
    if(!b) { return -1.0; }
    return out;
}

bool nearZero(double a) {
    if(a <= 0.000001 && a >= -0.000001) { return true; }
    return false;
}


class gen_map {
    public:
        vector<double> phys_pos;
        vector<double> gen_pos;
        gen_map(string filename, int idx, bool squish, bool cm) {
            ifstream gm;
            gm.open(filename.c_str());
            string line;
            double a,b;
            double rate = 1;
            if (cm) { rate = .01; }
            getline(gm,line);
            while(gm) {
                getline(gm,line);
                b = getGenCol(line,idx,a);
                b *= rate;
                if(gm.eof()) { break; }
                if (b <= -0.5 ) {
                    cerr << "Error with genetic map\n";
                    cerr << a << " " << b << " " << line << endl;
                    exit(1);
                }
                if(!squish || gen_pos.size() == 0 || (phys_pos.size() > 0 && gen_pos[gen_pos.size()-1] != b)) {
                    phys_pos.push_back(a);
                    gen_pos.push_back(b);
                }
                //if(gm.eof()) { break; }
            }
            cerr << phys_pos[0] << "\t" << gen_pos[0] << "\t" << phys_pos[phys_pos.size()-1] << "\t" << gen_pos[gen_pos.size()-1] << endl;
            cerr << phys_pos.size() << endl;
        }
        gen_map() { }
        double getGenPos(int position) {
            double dp = (double)position;
            double gr;
            if(dp < phys_pos[0]) {
                gr = (gen_pos[1]-gen_pos[0])/(phys_pos[1]-phys_pos[0]);
                return gr*(dp-phys_pos[0]);
            }
            int i = 1;
            while (i < phys_pos.size()-1 && dp > phys_pos[i]) { i++; }
            gr = (gen_pos[i]-gen_pos[i-1])/(phys_pos[i]-phys_pos[i-1]);
            return gen_pos[i-1]+gr*(dp-phys_pos[i-1]);
        }

};

vector<char> getAlleles(vector<string> vcfline) {
    vector<char> out;
    for(int i = 9; i < vcfline.size(); i++) {
        if (vcfline[i][0] != '0') { 
            out.push_back('1');
        } else { out.push_back('0'); }
        if (vcfline[i][2] != '0') {
            out.push_back('1');
        } else { out.push_back('0'); }
    }
    return out;
}

void printAD(vector<int> &a, vector<int> &d) {
    for (int i = 0; i < a.size(); i++) {
        cerr << d[i] << "\t";
    }
    cerr << endl;
    for(int i = 0; i < a.size(); i++) {
        cerr << a[i] << "\t";
    }
    cerr << endl;
}

int getAc(vector<char> alleles) {
    int out = 0;
    for(int i = 0; i < alleles.size(); i++) {
        if(alleles[i] == '1') {
            out += 1;
        }
    }
    return out;
}

snp_data parseVcfLine(string line) {
    vector<string> la = split(line,"\t");
    vector<char> alleles = getAlleles(la);
    int ac = getAc(alleles);
    int pos = atoi(la[1].c_str());
    snp_data out(ac,pos,alleles);
    return out;

}

void updateVectors(vector<int> & a, vector<int> & d, vector<char> alleles, int current_snp) {
    int p = current_snp+1;
    int q = current_snp+1;
    vector<int> a_new;
    vector<int> d_new;
    vector<int> b;
    vector<int> e;
    for(int i = 0; i < a.size(); i++) {
        p = max(d[i],p);
        q = max(d[i],q);
        if (alleles[a[i]] == '0') {
            a_new.push_back(a[i]);
            d_new.push_back(p);
            p = 0;
        } else {
            b.push_back(a[i]);
            e.push_back(q);
            q=0;
        }
    }
    a_new.insert(a_new.end(),b.begin(),b.end());
    d_new.insert(d_new.end(),e.begin(),e.end());
    for(int i = 0; i < a.size(); i++) {
        a[i] = a_new[i];
        d[i] = d_new[i];
    }
    return;
}

vector<string> msh(vector<int> & a, vector<int> &d, vector<int> &pos_list, int pos, int sample_count) {
    if (a.size() == 0) { 
        vector<string> out(sample_count,"0*");
        return out;
    }
    int l = a.size();
    vector<string> y_msh(l,"0");
    vector<string> ordered_msh(l,"0");
    int c_idx;
    for (int i = 0; i < l; i++) {
        if (i == l-1) { c_idx = d[l-1]; }
        else if (i == 0) { c_idx = d[1]; }
        else { c_idx = min(d[i],d[i+1]); }

        if (c_idx == 0) {
            int p = abs(pos-pos_list[0]);
            y_msh[i] = to_string(p)+'*';
        } else {
            int p = abs(pos - pos_list[c_idx-1]);
            y_msh[i] = to_string(p);
        }
    }
    for(int a_i = 0; a_i < a.size(); a_i++) {
        int a_v = a[a_i];
        ordered_msh[a_v] = y_msh[a_i];
    }
    return ordered_msh;

}

vector<int> getAPrimeVector(vector<int> &a){
    vector<int> aprime(a.size(),0);
    for(int i = 0; i < a.size();i++) {
        int a_v = a[i];
        aprime[a_v] = i;
    }
    return aprime;
}

/*int getDSingle(vector<int> &a, vector<int> &d, int idx, int sample_count) {
    if (a.size() == 0) { return -1; }
    //int *aidx;
    auto aidxit = find(a.begin(),a.end(),idx);
    auto aidx = distance(a.begin(),aidxit);
    if (aidx == sample_count - 1) { return d[sample_count-1]; }
    if (aidx == 0) { return d[1]; }
    return min(d[aidx],d[aidx+1]);
}*/

int findElementInVector(vector<int> &v, int i) {
    for (int j = 0; j < v.size(); j++) {
        if (v[j] == i) {
            return j;
        }
    }
    return -1;
}

int getDKton(vector<int> &a, vector<int> &aprime, vector<int> &d, int idx, int sample_count, vector<int> &idx_list) {
    if(aprime.size() == 0) { return -1; }
    int d_low = 0;
    int d_hi = 0;
    int a_idx = aprime[idx];
    bool hi_valid = false;
    bool low_valid = false;
    int i;
    if (a_idx != 0) {
        d_low = d[a_idx];
        low_valid = true;
    }
    if (a_idx != aprime.size()-1) {
        d_hi = d[a_idx+1];
        hi_valid = true;
    }
    i = a_idx-1;
    if (a_idx > 0 && findElementInVector(idx_list,a[a_idx-1]) != -1) {
        low_valid = false;
        while (i >= 1) {
            d_low = max(d_low,d[i]);
            if (findElementInVector(idx_list,a[i-1]) == -1) {
                low_valid = true;
                break;
            }
            i--;
        }
    }
    i = a_idx+2;
    if(a_idx < sample_count-1 && findElementInVector(idx_list,a[a_idx+1]) != -1) {
        hi_valid = false;
        while (i < sample_count) {
            d_hi = max(d_hi,d[i]);
            if (findElementInVector(idx_list,a[i]) == -1) {
                hi_valid = true;
                break;
            }
            i++;
        }
    }
    if (!hi_valid) { return d_low; }
    if (!low_valid) { return d_hi; }
    return min(d_hi,d_low);
}

int getDSingle(vector<int> &aprime, vector<int> &d, int idx, int sample_count) {
    if(aprime.size() == 0) { return -1; }
    int aidx = aprime[idx];
    if (aidx == sample_count - 1) { return d[sample_count-1]; }
    if (aidx == 0) { return d[1]; }
    return min(d[aidx],d[aidx+1]);
}

string getMshString(vector<int> &a, vector<int> &d, vector<int> &out_range,vector<int> &pos_list, vector<double> &gen_list, vector<int> &idx_list, int pos, double gpos, int sample_count, bool use_genmap, int round_sigfig) {
    stringstream outs;
    if (round_sigfig != -1) { outs << setprecision(round_sigfig); }
    outs << pos;

    if (use_genmap) { outs << "\t" << gpos; }
    vector<int> aprime = getAPrimeVector(a);
    for(int i = 0; i < out_range.size(); i++) {
        int d_idx;
        if (pos_list.size() == 0) {
            d_idx = -1;
        } else if (idx_list.size() == 0) {
            d_idx = getDSingle(aprime,d,out_range[i],sample_count);
        } else {
            d_idx = getDKton(a,aprime,d,out_range[i],sample_count,idx_list);
        }
        //auto aidxit = find(a.begin(),a.end(),i);
        //auto aidx = distance(a.begin(),aidxit);
        outs << "\t";// << d_idx << "\t" << aprime[i] << "\t";
        if (d_idx == -1) { outs << "0*"; }
        else if (d_idx == 0) {
            outs << abs(pos-pos_list[d_idx]) << "*";
        } else {
            outs << abs(pos-pos_list[d_idx-1]);
        }
        if (use_genmap) {
            if(d_idx == -1) { outs << ":0.0*"; }
            else if (d_idx == 0) { 
                outs << ":" << abs(gpos-gen_list[0]) << "*";
            } else {
                outs << ":" << abs(gpos-gen_list[d_idx-1]);
            }
        }

    }
    outs << "\n";
    return outs.str();
}

vector<int> getKIdxList(snp_data &sd) {
    //char target_allele = '0';
    vector<int> out_idx;
    /*if(sd.allele_count == k) {
        target_allele = '1';
    }*/
    for(int i = 0; i < sd.alleles.size(); i++) {
        if (sd.alleles[i] == '1') {
            out_idx.push_back(i);
        }
    }
    return out_idx;
}

int getSingletonIdx(snp_data sd) {
    int o;
    if (sd.allele_count == 1) {
        o = *find(sd.alleles.begin(),sd.alleles.end(),'1');
    } else {
        o = *find(sd.alleles.begin(),sd.alleles.end(),'0');
    }
    return (o/2)*2;
}

vector<int> getSingletonIdxList(snp_data &sd) {
    vector<char>::iterator op;
    int o;
    if (sd.allele_count == 1) {
        op = find(sd.alleles.begin(),sd.alleles.end(),'1');
    } else {
        op = find(sd.alleles.begin(),sd.alleles.end(),'0');
    }
    o = distance(sd.alleles.begin(),op);
    o /= 2;
    o *= 2;
    vector<int> out;
    out.push_back(o);
    out.push_back(o+1);
    return out;
}

vector<int> getFullIdxList(int sample_count) {
    vector<int> out;
    for(int i = 0; i < sample_count; i++) {
        out.push_back(i);
    }
    return out;
}

vector<int> a_start(int sc) {
    vector<int> out;
    for (int i = 0; i < sc; i++) {
        out.push_back(i);
    }
    return out;
}

vector<int> d_start(int sc) {
    vector<int> out(sc,0);
    return out;
}

int main(int argc, char ** argv) {
    setbuf(stdout,NULL);
    string filename = "-";// = argv[1];
    bool include_singletons = false;
    bool singleton_mode = false;
    string outfilename = "default.out";
    string mapname = "";
    int map_colidx = 0;
    int round_sigfig = -1;
    bool squish_map = true;
    bool map_cm = true;
    bool use_genmap = false;
    bool k_mode = false;
    bool k_all = false;
    int k_low = -1;
    int k_hi = -1;
    for(int i = 1; i < argc; i++) {
        string arg = argv[i];
        if(arg == "--vcf") { filename = argv[++i]; }
        else if(arg == "--include-singletons") { include_singletons = true; }
        else if(arg == "--singleton") { singleton_mode = true; }
        else if(arg == "--out") { outfilename = argv[++i]; }
        else if(arg == "--map") { use_genmap = true; mapname = argv[++i]; }
        else if(arg == "--gen-idx") { map_colidx = atoi(argv[++i]); }
        else if(arg == "--nosquish") { squish_map = false; }
        else if(arg == "--round") { round_sigfig = atoi(argv[++i]); }
        else if(arg == "--k-all") { k_mode = true; k_all = true; }
        else if(arg == "--k-range") { k_mode = true; k_low = atoi(argv[++i]); 
                                      k_hi = atoi(argv[++i]);
        }
        else { cerr << "Argument " << arg << " not recognized\n"; exit(1); }
    }
    istream *ins;
    ifstream infs;
    if (filename.compare("-") == 0) {
        ins = &cin;
    } else {
        infs.open(filename.c_str());
        ins = &infs;
    }
    gen_map genetic_map;
    if (mapname.compare("") != 0) {
        genetic_map = gen_map(mapname,map_colidx,squish_map,map_cm);
    }
    if (round_sigfig != -1) {
        cout << setprecision(round_sigfig);
    }
    string line;
    int snp_count = 0;
    int sample_count = 0;
    vector<int> a;
    vector<int> d;
    vector<int> pos_list;
    vector<double> gen_list;
    ofstream outf;
    outf.open(outfilename);

    while (getline(*ins,line)) {
        if (line.size() == 0) { break; }
        if (line[0] == '#') { continue; }
        vector<string> la = split(line,"\t");
        if (la[3].size() != 1 || la[4].size() != 1) { continue; }
        if (la[7].find("CPG_TAG") != string::npos) { continue; }
        snp_data sd = parseVcfLine(line);
        bool is_singleton = false;
        int noninf_pos = -1;
        int k_pos = -1;
        double noninf_gen, k_gen;
        vector<int> out_range;
        vector<int> k_idxlist;
        if (sd.allele_count == 0 || sd.allele_count == sd.alleles.size()) { continue; }
        if (sd.allele_count == 1 || sd.allele_count == sd.alleles.size()-1 ) {
            if (!singleton_mode && !k_mode) {
                if (!include_singletons) {
                    continue;
                }
            } else if (!k_mode) {

                noninf_pos = stoi(la[1]);
                if(use_genmap) { noninf_gen = genetic_map.getGenPos(noninf_pos); }
                //singleton_idx = getSingletonIdx(sd);
            }
        }
        if (k_mode && (k_all || (sd.allele_count >= k_low && sd.allele_count <= k_hi))) {
            k_idxlist = getKIdxList(sd);
            k_pos = atoi(la[1].c_str());
            if (use_genmap) {
                k_gen = genetic_map.getGenPos(k_pos);
            }
        }
        if (snp_count == 0) {
            sample_count = sd.alleles.size();
            a = a_start(sample_count);
            d = d_start(sample_count);
        }
        if(singleton_mode && (noninf_pos != -1)) {
            out_range = getSingletonIdxList(sd);
            string out_string = getMshString(a,d,out_range,pos_list,gen_list,k_idxlist,noninf_pos,noninf_gen,sample_count,use_genmap,round_sigfig);
            outf << out_string;
        }
        if(k_mode) {
            string out_string = getMshString(a,d,k_idxlist,pos_list,gen_list,k_idxlist,k_pos,k_gen,sample_count,use_genmap,round_sigfig);
            outf << out_string;
        }
        if (!singleton_mode || noninf_pos == -1 || include_singletons) {
            int pos_add = atoi(la[1].c_str());
            pos_list.push_back(pos_add);
            if(use_genmap) { gen_list.push_back(genetic_map.getGenPos(pos_add)); }
            updateVectors(a,d,sd.alleles,snp_count);
            snp_count += 1;
        }
        if (!singleton_mode && !k_mode) {
            out_range = getFullIdxList(sample_count);
            double gpos = (gen_list.size() == 0 ? 0 : gen_list[gen_list.size()-1]);
            string out_string = getMshString(a,d,out_range,pos_list,gen_list,k_idxlist,pos_list[pos_list.size()-1],gpos,sample_count,use_genmap,round_sigfig);
            outf << out_string;
        }
    }

    return 0;
}