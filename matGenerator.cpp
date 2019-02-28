#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <array>

using namespace std;

void readMatrix(const string &filename, vector<vector<int>> &M, bool binary = false) {
    ifstream inFile (filename);
    string line;
    char del = '\t';

    if(!inFile)
    {
        cerr << "Error, file does not exist. " << endl;
        exit(EXIT_FAILURE);
    }

    while(getline(inFile, line)){
        istringstream tokenStream(line);
        string token;
        vector<int> tokens;
        while(getline(tokenStream, token, del)){
            if (binary)
                tokens.push_back(bool(stod(token)));
            else
                tokens.push_back(stoi(token));
        }
        M.push_back(tokens);
    }

    inFile.close();
}

void readMatrix(const string &filename, vector<vector<double>> &M) {
    ifstream inFile (filename);
    string line;
    char del = '\t';

    if(!inFile)
    {
        cerr << "Error, file does not exist. " << endl;
        exit(EXIT_FAILURE);
    }

    while(getline(inFile, line)){
        istringstream tokenStream(line);
        string token;
        vector<double> tokens;
        while(getline(tokenStream, token, del)){
            tokens.push_back(stod(token));
        }
        M.push_back(tokens);
    }

    inFile.close();
}

// given vaf matrix
// p_i,j is primarily indexed by j
void writeMatrix(const string &filename, vector<vector<int>> &M, vector<vector<double>> &V, const int &K, const int &T) {
    // define other parameters
    auto m = static_cast<int>(M.size());
    auto n = static_cast<int>(M[0].size());

    // construct constraint matrix
    ofstream outFile (filename);
    if(!outFile)
    {
        cerr << "Error, cannot write into file. " << endl;
        exit(EXIT_FAILURE);
    }

    // dimension and element number
    int nrow = n + K + m + T + m*n*K + m*K + T*K + m*T*K*(K-1) + m*n;
    int ncol = n*K+m*n*K+m*T+T*K*K;
    int nele = n*K + K*n + m*T + T*m + m*n*K*2 + m*K*n + 2*T*K*(K-1) + m*T*K*(K-1)*(2*n+2) + m*n*K;
    int nzero = 0;
    for (int i = 0; i < V.size(); i++)
        for (int j = 0; j < V[0].size(); j++)
            if (V[i][j] == 0)
                nzero++;
    nele -= nzero*2*T*K*(K-1);

    outFile << nrow << " " << ncol << " " << nele << endl;

    int basej;

    // pathway_uniqueness "each column is assigned to exactly one set"
    basej = 0;
    for (int j = 1; j <= n; j++) {
        for (int k = 1; k <= K; k++)
            outFile << basej+ (j-1)*K + k << " " << 1 << " ";
        outFile << endl;
    }

    // pathway_nonempty "for each set Pk, at least one column is assigned to it"
    basej = 0;
    for (int k = 1; k <= K; k++) {
        for (int j = 1; j <= n; j++)
            outFile << basej+ (j-1)*K + k << " " << 1 << " ";
        outFile << endl;
    }


    // subtype_uniqueness "each row is assigned to exact one tree"
    basej = n*K+m*n*K;
    for (int i = 1; i <= m; i++) {
        for (int t = 1; t <= T; t++)
            outFile << basej+ (i-1)*T + t << " " << 1 << " ";
        outFile << endl;
    }

    // subtype_nonempty "for each subtype, at least one row is assigned to it"
    basej = n*K+m*n*K;
    for (int t = 1; t <= T; t++) {
        for (int i = 1; i <= m; i++)
            outFile << basej+ (i-1)*T + t << " " << 1 << " ";
        outFile << endl;
    }

    // necessary condition of mutation "gene can be mutated in pathway k only if it belongs to pathway k"
    basej = 0;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            for (int k = 1; k <= K; k++)
                outFile << basej+ (j-1)*K + k << " " << -1 << " " << basej+n*K + (i-1)*n*K + (j-1)*K + k << " " << 1 << endl;
        }
    }

    // mutual exclusivity "each patient has at most one mutation in one pathway"
    basej = n*K;
    for (int i = 1; i <= m; i++) {
        for (int k = 1; k <= K; k++) {
            for (int j = 1; j <= n; j++)
                outFile << basej+ (i-1)*n*K + (j-1)*K + k << " " << 1 << " ";
            outFile << endl;
        }
    }

    // linear progression constraint
    basej = n*K + m*n*K + m*T;
    for (int t = 1; t <= T; t++) {
        for (int k = 1; k <= K; k++) {
            for (int l = 1; l <= K; l++)
                if (l != k){
                    outFile << basej+ (t-1)*K*K + (k-1)*K + l << " " << 1 << " ";
                    outFile << basej+ (t-1)*K*K + (l-1)*K + k << " " << 1 << " ";
                }
            outFile << endl;
        }
    }

    // progression constraint "each patient satisfies the progression model of its SPM"
    basej = n*K;
    for (int i = 1; i <= m; i++)
        for (int t = 1; t <= T; t++)
            for (int k = 1; k <= K; k++)
                for (int l = 1; l <= K; l++)
                    if (l != k) {
                        outFile << basej+m*n*K+ (i-1)*T + t << " " << -1 << " ";
                        outFile << basej+m*n*K+m*T+ (t-1)*K*K + (k-1)*K + l << " " << -1 << " ";
                        for (int j = 1; j <= n; j++) {
                            if (V[i-1][j-1] > 0) {
                                outFile << basej+ (i-1)*n*K + (j-1)*K + k << " " << V[i-1][j-1] << " ";
                                outFile << basej+ (i-1)*n*K + (j-1)*K + l << " " << -V[i-1][j-1] << " ";
                            }
                        }
                        outFile << endl;
                    }

    // no false negative assumption "no real mutation is missed in observation"
    basej = n*K;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            for (int k = 1; k <= K; k++)
                outFile << basej+ (i-1)*n*K + (j-1)*K + k << " " << 1 << " ";
            outFile << endl;
        }
    }

    outFile.close();
}

// p_i,j is primarily indexed by j
void writePLPM(const string &filename, vector<vector<int>> &M, const int &K) {
    // define other parameters
    auto m = static_cast<int>(M.size());
    auto n = static_cast<int>(M[0].size());

    // construct constraint matrix
    ofstream outFile (filename);
    if(!outFile)
    {
        cerr << "Error, cannot write into file. " << endl;
        exit(EXIT_FAILURE);
    }

    // dimension and element number
    int nrow = n + K + m*(K-1) + m*K;
    int ncol = n*K + m*K + m*K;
    int nele = n*K + n*K + 2*m*(K-1) + m*K*(n+2);
    int nzero = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (M[i][j] == 0)
                nzero++;
    nele -= nzero*K;

    outFile << nrow << " " << ncol << " " << nele << endl;

    int basej;

    // constraint "each column is assigned to exactly one set"
    basej = 0;
    for (int j = 1; j <= n; j++) {
        for (int k = 1; k <= K; k++)
            outFile << basej+ (j-1)*K + k << " " << 1 << " ";
        outFile << endl;
    }

    // constraint "for each set Pk, at least one column is assigned to it"
    basej = 0;
    for (int k = 1; k <= K; k++) {
        for (int j = 1; j <= n; j++)
            outFile << basej+ (j-1)*K + k << " " << 1 << " ";
        outFile << endl;
    }


    // constraint "for each sample the progression model is satisfied"
    basej = n*K;
    for (int i = 1; i <= m; i++) {
        for (int k = 1; k < K; k++) {
            outFile << basej + (i - 1) * K + k << " " << 1 << " " << basej + (i - 1) * K + k + 1 << " " << -1 << " ";
            outFile << endl;
        }
    }

    // constraint "for each row ri, the set Pk is considered mutated if it has a 1 in ri
    // or if one of its entries in row ri is flipped to make it mutated"
    basej = 0;
    for (int i = 1; i <= m; i++) {
        for (int k = 1; k <= K; k++) {
            for (int j = 1; j <= n; j++)
                if (M[i-1][j-1] > 0)
                    outFile << basej+ (j-1)*K + k << " " << M[i-1][j-1] << " ";

            outFile << basej+n*K+ (i-1)*K + k << " " << -1 << " ";
            outFile << basej+n*K+m*K+ (i-1)*K + k << " " << 1 << " ";
            outFile << endl;
        }
    }

    outFile.close();
}

//Generate sparse matrix for SPM
void SPM_Sparse() {
    //user defined parameters
    //--------------------------
    // K is the expected size of pathways
    // T is the expected number of subtypes
    int K = 4;
    int T = 2;
    // file names of input and output
    string in = "GBM_SPM_CCF_Matrix.txt";
    string out = "anc_K" + to_string(K) + "_T" + to_string(T) + "_" + in;
    // directories of input and output
    string addrin = "../../data/COMB/";
    string addrout = "../../results/SPMInput/GBM_Ancestry/";
    //--------------------------

    // vaf matrix V and binarized matrix M
    vector<vector<int>> M;
    vector<vector<double>> V;

    readMatrix(addrin+in, V);
    readMatrix(addrin+in, M, true);

    // generate matrix file
    writeMatrix(addrout+out, M, V, K, T);
}

//Generate sparse matrix for PLPM
void PLPM_Sparse() {
    int K = 5;
    // file names of input and output
    string in = "WOOD_SPM_Binary_Matrix.txt";
    string out = "PLPM_K" + to_string(K) + "_" + in;
    // directories of input and output
    string addrin = "../../data/COMB/";
    string addrout = "../../results/PLPMInput/PLPM_WOOD/";

    // vaf matrix V and binarized matrix M
    vector<vector<int>> M;

    readMatrix(addrin+in, M, true);

    // generate matrix file
    writePLPM(addrout+out, M, K);
}


int main() {

    SPM_Sparse();
    return 0;
}
