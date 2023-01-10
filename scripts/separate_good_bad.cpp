#include <iostream>
#include <map>
#include <fstream>
#include <string> 
#include <sstream>
#include <queue> 
#include <thread>
#include <unistd.h>
#include <mutex>
#include <ctype.h>
using namespace std;

class SeqCarrier {
    public:
        string read1;
        string read2;
        string reference;
        SeqCarrier(){}
        SeqCarrier(string r1, string r2, string ref){
            read1 = r1;
            read2 = r2;
            reference = ref;
        }
};

mutex mtx;
queue<SeqCarrier> seq_q;
bool running = true;

bool all_in_k(vector<bool>& b, int k, int ml, int st, int nd){
    bool out = true;
    for(int i = st; i < nd; i++){
        out &= b[i];
    }
    return out;
}
void separate_good_bad(string path, string upath, int k){
    SeqCarrier seqs;
    ofstream out(path);
    ofstream unknown(upath);
    int ml, st, nd;
    bool b1, b2;
    vector<bool> r1_eq_ref(150);
    vector<bool> r2_eq_ref(150);
    vector<bool> r1_eq_r2(150);
    while(running){
        if(seq_q.size() > 0){
            mtx.lock();
            seqs = seq_q.front();
            seq_q.pop();
            mtx.unlock();
            ml = min(seqs.read1.length(),min(seqs.read2.length(),seqs.reference.length()));
            //cout << seqs.reference  << endl;
            for(int i=0; i < ml; i++){
                r1_eq_ref[i] = seqs.read1[i] == toupper(seqs.reference[i]);
                r2_eq_ref[i] = seqs.read2[i] == toupper(seqs.reference[i]);
                r1_eq_r2[i] = seqs.read1[i] == seqs.read2[i];
            }
            for(int jj=0; jj < ml; jj++){
             //   cout << r1_eq_ref[i] << " " << r2_eq_ref[i] << r1_eq_r2[i] << endl;
                if(r1_eq_ref[jj] && ! r2_eq_ref[jj]){
                    st = max(0, jj-(k/2));
                    b1 = all_in_k(r1_eq_ref, k, ml, st, nd);
                    if(b1 && (st >= k/2)){
                        out << seqs.read1.substr(st, k) << "\t" << seqs.read2.substr(st, k) << endl; //<< "\t" << seqs.read1  
                        //    << "\n%" << seqs.read2.substr(st, k) << "\t" << seqs.read2
                        //    << "\n@" << seqs.reference.substr(st,k) << "\t" << seqs.reference << endl;
                    }
                }
                if( ! r1_eq_ref[jj] && r2_eq_ref[jj]){
                    st = max(0, jj-(k/2));
                    b1 = all_in_k(r2_eq_ref, k, ml, st, nd);
                    if(b1 && (st >= k/2)){
                        out << seqs.read2.substr(st, k) << "\t"  << seqs.read1.substr(st, k) << endl;
                    }    
                }
                if ( ! r1_eq_ref[jj] && r1_eq_r2[jj] ){
                    st = max(0, jj-(k/2));
                    b1 = all_in_k(r1_eq_r2, k, ml, st, nd);
                    if(b1 && (st >= k/2)){
                        unknown << seqs.reference.substr(st, k) << "\t" << seqs.read1.substr(st, k) << endl;
                    }    
                }
            }
        } else {
            usleep(5);//wait for data
        }
    }
    out.close();
    unknown.close();
}

void read_overlaps(string path, map<string, string> ref_fasta){
    cout << "Read overlaps" << endl;
    map<string, SeqCarrier> output;
    ifstream input;
    input.open(path);
    if(!input.is_open()){
        cerr << "Failed to open " << path << endl;
    }
    string line;
    string piece;
    int pos;
    stringstream linebuffer;
    vector<string> pieces(10);
    SeqCarrier s;
    string ref;
    bool flag;
    while(getline(input, line)){
        linebuffer.str(line);
        linebuffer.clear();
        pieces.clear();
        while(getline(linebuffer, piece, '\t')){
            pieces.push_back(piece);
        }
        if(pieces.size() == 9){
            pos = stoi(pieces[6])-1;
            //cout << pieces[3] << endl;
            flag = false;
            char cx;
            // if there's anything other than a match/mismatch, continue
            for(char cx : pieces[1]){
                flag |= (cx == 'D' | cx == 'I' | cx == 'S');
            } 
            for(char cx : pieces[2]){
                flag |= (cx == 'D' | cx == 'I' | cx == 'S');
            }
            if(flag) continue;
            if(auto result = ref_fasta.find(pieces[3]); result != ref_fasta.end()){
                if(pos > result->second.length()){
                     cerr << pieces[3] << " " << result->second.length() << " " << pos << endl;
                     continue;
                }
                ref = result->second.substr(pos, pieces[7].length());
                //cout << pos << " " << pieces[7].length() << " " << ref << endl;
                s = SeqCarrier(pieces[7], pieces[8], ref);
                mtx.lock();
                seq_q.push(s);
                mtx.unlock();
            } else {
                continue;
            }
        }
    }
    return;
} 


map<string, string> parse_fasta(string path){
    map<string, string> fasta;
    string line;
    ifstream input;
    input.open(path);
    if(!input.is_open()) {
        cerr << "Failed to open " << path << endl;
        return fasta;
    }
    stringstream buffer;
    string key;
    string sequence;
    while(getline(input, line)){
        if(line[0] == '>'){
            sequence = buffer.str();
            if (!sequence.empty()){
                fasta.insert(pair<string, string>(key, sequence));
                buffer.str("");
            }
            key = line.substr(1);
        } else{
            buffer << line;
        }
    }
    sequence = buffer.str();
    if(!sequence.empty()){
        fasta.insert(pair<string, string>(key, sequence));
    }
    return fasta;

}

int main(int argc, char** argv){
    string path = argv[1]; //"chr1.fa";
    string overlaps = argv[2];//"new_new_overlaps.txt";
    string out = argv[3]; //"subs_good.txt";
    string unknown = argv[4]; //"subs_unknown.txt";
    int k = stoi(argv[5]);  //"17");
    map<string, string> fa = parse_fasta(path);
    for(auto p : fa) cout << p.first << endl;
    thread reader(read_overlaps,overlaps, fa);
    thread separator(separate_good_bad, out, unknown, k);
    reader.join();
    running = false;
    separator.join();
    return 0;
}
