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
            //cerr << position << "\t" << allele_count << endl;
        }
        snp_data() {
            allele_count = -1;
            position = -1;
        }
};

double getGenCol(string line, int col, double & phys) {
    stringstream s(line);
    s >> phys;
    double j = 0, out = 0;
    for(int i = 0; i < col; i++) {
        s >> j;
    }
    bool b = bool(s >> out);
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
            //getline(gm,line);
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
            //cerr << phys_pos[0] << "\t" << gen_pos[0] << "\t" << phys_pos[phys_pos.size()-1] << "\t" << gen_pos[gen_pos.size()-1] << endl;
            //cerr << phys_pos.size() << endl;
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

//snp_data parseVcfLine(string & line) {
//    vector<string> la = split(line,"\t");
snp_data parseVcfLine(vector<string> &la) {
    vector<char> alleles = getAlleles(la);
    int ac = getAc(alleles);
    int pos = atoi(la[1].c_str());
    snp_data out(ac,pos,alleles);
    return out;

}

snp_data parseVcfLine2(string & line, int presize, vector<int> &sub_list, bool sub_flag) {
    stringstream s(line);
    string junk, info, ref, alt, chrom;
    int phys;
    vector<char> alleles(presize,'0');
    s >> chrom >> phys >> junk >> ref >> alt >> junk >> junk >> info >> junk;
    if (ref.size() != 1 || alt.size() != 1 || info.find("CPG_TAG") != string::npos) {
        snp_data out;
        return out;
    }
    string hap;
    int ac = 0;
    int sub_idx = 0;
    int sub_total = sub_list.size();
    int hap_count = 0;
    while (s >> hap) {
        /*if (hap[0] != '0') {
            ac++;
        }
        if (hap[2] != '0') {
            ac++;
        }
        alleles.push_back(hap[0]);
        alleles.push_back(hap[2]);*/
        int geno_idx = 0;
        string delim_str = "/|";
        while (geno_idx == 0 || delim_str.find(hap[geno_idx-1]) != string::npos) {
            if (sub_flag) {
                if (sub_list[sub_idx] == hap_count) {
                    sub_idx++;
                } else {
                    geno_idx+=2;
                    hap_count++;
                    continue;
                }
            }
            if (hap[geno_idx] != '0') {
                ac++;
            }
            alleles.push_back(hap[geno_idx]);
            geno_idx+=2;
            hap_count++;
            if (sub_total != 0 && sub_total == sub_idx) { break; }
        }

    }
    snp_data out(ac,phys,alleles);
    return out;
}

void updateVectors(vector<int> & a, vector<int> & d, vector<char> alleles, int current_snp) {
    int p = current_snp+1;
    int q = current_snp+1;
    int size = a.size();
    vector<int> a_new(size,0);
    vector<int> d_new(size,0);
    vector<int> b(size,0);
    vector<int> e(size,0);
    int ref_i = 0;
    int alt_i = 0;
    for(int i = 0; i < size; i++) {
        p = max(d[i],p);
        q = max(d[i],q);
        if (alleles[a[i]] == '0') {
            //a_new.push_back(a[i]);
            //d_new.push_back(p);
            a_new[ref_i] = a[i];
            d_new[ref_i] = p;
            ref_i++;
            p = 0;
        } else {
            //b.push_back(a[i]);
            //e.push_back(q);
            b[alt_i] = a[i];
            e[alt_i] = q;
            alt_i++;
            q=0;
        }
    }
    /*a_new.insert(a_new.end(),b.begin(),b.end());
    d_new.insert(d_new.end(),e.begin(),e.end());
    for(int i = 0; i < a.size(); i++) {
        a[i] = a_new[i];
        d[i] = d_new[i];
    }*/
    int alt_size = ref_i;
    for(int i = 0; i < alt_size; i++) {
        a[i] = a_new[i];
        d[i] = d_new[i];
    }
    for(int i = alt_size; i < size; i++) {
        a[i] = b[i-alt_size];
        d[i] = e[i-alt_size];
    }
    return;
}

vector<string> msh(vector<int> & a, vector<int> &d, vector<int> &pos_list, int pos, int sample_count, bool print_indiv) {
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

int getASingle(vector<int> &aprime, vector<int> &a, vector<int> &d, int idx, int sample_count) {
    if(aprime.size() == 0) { return -1; }
    int aidx = aprime[idx];
    if (aidx == sample_count - 1) { return a[sample_count-2]; }
    if (aidx == 0) { return a[1]; }
    if (d[aidx] < d[aidx+1]) { return a[aidx-1]; }
    else { return a[aidx+1]; }
}

string getMshString(vector<int> &a, vector<int> &d, vector<int> &out_range,vector<int> &pos_list, vector<double> &gen_list, vector<int> &idx_list, int pos, double gpos, int sample_count, bool use_genmap, int round_sigfig, bool print_chrom) {
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
        if (print_chrom) {
            int this_indiv = out_range[i];
            int msh_indiv = getASingle(aprime,a,d,out_range[i],sample_count);
            if (d_idx <= 0) {
                outs << ":-1:-1";
            } else {
                outs << ":" << this_indiv << ":" << msh_indiv;
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

bool positionCondition(vector<int> &outpos_list, int outpos_idx, int pos, bool revpos_flag) {
    if (outpos_idx == outpos_list.size()) { return false; }
    if (revpos_flag && outpos_list[outpos_idx] >= pos) { return true; }
    if (!revpos_flag && outpos_list[outpos_idx] < pos) { return true; }
    return false; 
}

int main(int argc, char ** argv) {
    if (argc == 1 || argv[1] == "-h" || argv[1] == "--help") {
        cerr << "usage: ./msh_vcf [-h] [--vcf VCFNAME] [--gen GENNAME] [--sub SUBNAME]\n"
        "               [--out OUTNAME] [--gen-idx GENIDX] [--nosquish]\n"
        "               [--round ROUND]\n"
        "               [--k-all | --k-range K_START K_END | --singleton]\n"
        "               [--positions POSNAME]\n"
        "               [--exclude-singletons] [--revpos]\n"
        "\n"
        "Arguments:\n"
        "-h, --help            show this help message and exit\n"
        "--vcf VCFNAME         Input VCF filename, either uncompressed or gzipped\n"
        "--gen GENNAME         Name of map file\n"
        "--sub SUBNAME         Name of file with subsample indices\n"
        "--out OUTNAME         Name of output msh file\n"
        "--gen-idx GENIDX      Use [i+2]th column of genetic map file for genetic\n"
        "                      distances\n"
        "--nosquish            Read all rows in map file, not just rows with\n"
        "                      differing genetic distances\n"
        "--round ROUND         Round floats to this many significant figures\n"
        "--k-all               Return estimatesfor every variant\n"
        "--k-range K_START K_END\n"
        "                      Return estimates for every variant with allele count\n"
        "                      in given range\n"
        "--singleton           Singleton mode: Generate two msh values, one for each\n"
        "                      haplotype of an individual with a singleton\n"
        "--positions POSNAME   File with positions for alledges estimates\n"
        "--exclude-singletons  DO NOT Allow singletons to terminate MSH lengths\n"
        "--revpos              If using --positions, indicate this is for right\n"
        "                      length generation and input VCF is reversed\n";
        exit(1);

    }
    setbuf(stdout,NULL);
    string filename = "-";// = argv[1];
    bool include_singletons = true;
    bool singleton_mode = false;
    string outfilename = "default.out";
    string mapname = "";
    string posname = "";
    string subname = "";
    int map_colidx = 0;
    int round_sigfig = -1;
    bool squish_map = true;
    bool map_cm = true;
    bool use_genmap = false;
    bool pos_flag = false;
    bool k_mode = false;
    bool k_all = false;
    bool revpos = false;
    bool sub_flag = false;
    bool print_indiv = false;
    int k_low = -1;
    int k_hi = -1;
    for(int i = 1; i < argc; i++) {
        string arg = argv[i];
        if(arg == "--vcf") { filename = argv[++i]; }
        else if(arg == "--exclude-singletons") { include_singletons = false; }
        else if(arg == "--singleton") { singleton_mode = true; }
        else if(arg == "--out") { outfilename = argv[++i]; }
        else if(arg == "--gen") { use_genmap = true; mapname = argv[++i]; }
        else if(arg == "--gen-idx") { map_colidx = atoi(argv[++i]); }
        else if(arg == "--nosquish") { squish_map = false; }
        else if(arg == "--round") { round_sigfig = atoi(argv[++i]); }
        else if(arg == "--positions") { pos_flag = true; posname = argv[++i]; }
        else if(arg == "--revpos") { revpos = true; }
        else if(arg == "--sub") { sub_flag = true; subname = argv[++i]; }
        else if(arg == "--k-all") { k_mode = true; k_all = true; }
        else if(arg == "--k-range") { k_mode = true; k_low = atoi(argv[++i]); 
                                      k_hi = atoi(argv[++i]);
        }
        else if (arg == "--print-indiv") { print_indiv = true; }
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
    vector<int> outpos_list;
    ifstream posf;
    string line;
    int tpos;
    int outpos_idx;
    if (posname.compare("") != 0) {
        posf.open(posname.c_str());
        while(getline(posf,line)) {
            tpos = atoi(line.c_str());
            if (revpos) {
                outpos_list.insert(outpos_list.begin(),tpos);
            } else {
                outpos_list.push_back(tpos);
            }
        }
        outpos_idx = 0;
        pos_flag = true;
    }

    vector<int> sub_list;
    ifstream subf;
    if (subname.compare("") != 0) {
        subf.open(subname.c_str());
        while (getline(subf,line)) {
            tpos = atoi(line.c_str());
            sub_list.push_back(tpos);
        }
        sort(sub_list.begin(),sub_list.end());
        sub_flag = true;
    }

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
        /*vector<string> la = split(line,"\t");
        if (la[3].size() != 1 || la[4].size() != 1) { continue; }
        if (la[7].find("CPG_TAG") != string::npos) { continue; }
        snp_data sd = parseVcfLine(la);*/
        snp_data sd = parseVcfLine2(line,0,sub_list,sub_flag);
        if (sd.position == -1) { continue; }
        //cerr << sd.allele_count << "\t" << sd.position << endl;
        bool is_singleton = false;
        int noninf_pos = -1;
        int k_pos = -1;
        double noninf_gen, k_gen;
        vector<int> out_range;
        vector<int> k_idxlist;
        if (sample_count == 0) { sample_count = sd.alleles.size(); }
        if (sd.allele_count == 0 || sd.allele_count == sample_count) { continue; }
        if (sd.allele_count == 1 || sd.allele_count == sample_count-1 ) {
            if (!singleton_mode && !k_mode) {
                if (!include_singletons) {
                    continue;
                }
            } else if (!k_mode) {

                noninf_pos = sd.position;//stoi(la[1]);
                if(use_genmap) { noninf_gen = genetic_map.getGenPos(noninf_pos); }
                //singleton_idx = getSingletonIdx(sd);
            }
            if (k_all && sd.allele_count == 1 && !include_singletons) {
                continue;
            }
        }
        bool do_k = (k_all && sd.allele_count != 1);

        if (k_mode && (do_k || (sd.allele_count >= k_low && sd.allele_count <= k_hi))) {
            k_idxlist = getKIdxList(sd);
            k_pos = sd.position;//atoi(la[1].c_str());
            if (use_genmap) {
                k_gen = genetic_map.getGenPos(k_pos);
            }
        }
        if (snp_count == 0) {
            //sample_count = sd.alleles.size();
            a = a_start(sample_count);
            d = d_start(sample_count);
        }
        if (pos_flag) {
            while(positionCondition(outpos_list,outpos_idx,sd.position,revpos)) {
                int opos = outpos_list[outpos_idx];
                out_range = getFullIdxList(sample_count);
                double gpos = (gen_list.size() == 0 ? 0 : genetic_map.getGenPos(opos));
                string out_string = getMshString(a,d,out_range,pos_list,gen_list,k_idxlist,opos,gpos,sample_count,use_genmap,round_sigfig,print_indiv);
                outf << out_string;
                outpos_idx++;
            }
        }
        if(singleton_mode && (noninf_pos != -1)) {
            out_range = getSingletonIdxList(sd);
            string out_string = getMshString(a,d,out_range,pos_list,gen_list,k_idxlist,noninf_pos,noninf_gen,sample_count,use_genmap,round_sigfig,print_indiv);
            outf << out_string;
        }
        if(k_pos != -1) {
            string out_string = getMshString(a,d,k_idxlist,pos_list,gen_list,k_idxlist,k_pos,k_gen,sample_count,use_genmap,round_sigfig,print_indiv);
            outf << out_string;
        }
        if (!singleton_mode || noninf_pos == -1 || include_singletons) {
            int pos_add = sd.position;//atoi(la[1].c_str());
            pos_list.push_back(pos_add);
            if(use_genmap) { gen_list.push_back(genetic_map.getGenPos(pos_add)); }
            updateVectors(a,d,sd.alleles,snp_count);
            snp_count += 1;
        }
        if (!singleton_mode && !k_mode && !pos_flag) {
            out_range = getFullIdxList(sample_count);
            double gpos = (gen_list.size() == 0 ? 0 : gen_list[gen_list.size()-1]);
            string out_string = getMshString(a,d,out_range,pos_list,gen_list,k_idxlist,pos_list[pos_list.size()-1],gpos,sample_count,use_genmap,round_sigfig,print_indiv);
            outf << out_string;
        }
    }

    return 0;
}