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
    int nrow = n + K + m + T + m*n*K + m*K + 2*T*(K+1) + m*T*K*(K-1) + m*n;
    int ncol = n*K+m*n*K+m*T+T*(K+1)*(K+1);
    int nele = n*K + K*n + m*T + T*m + m*n*K*2 + m*K*n + 2*T*(K+1)*K + m*T*K*(K-1)*(2*n+2) + m*n*K;
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
        for (int k = 1; k <= K+1; k++) {
            for (int l = 1; l <= K+1; l++)
                if (l != k)
                    outFile << basej+ (t-1)*(K+1)*(K+1) + (k-1)*(K+1) + l << " " << 1 << " ";
            outFile << endl;

            for (int l = 1; l <= K+1; l++)
                if (l != k)
                    outFile << basej+ (t-1)*(K+1)*(K+1) + (l-1)*(K+1) + k << " " << 1 << " ";
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
                        outFile << basej+m*n*K+m*T+ (t-1)*(K+1)*(K+1) + (k-1)*(K+1) + l << " " << -1 << " ";
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

int main() {
    //-----------CONCISE_MATRIX-------
    for (int K = 3; K <= 7; K++)
        for (int T = 1; T <= 3; T++) {
            // file names of input and output
            string in = "TCGA_COAD_CCFs_spm.txt";
            string out = "ccs_K" + to_string(K) + "_T" + to_string(T) + "_" + in;
            // directories of input and output
            string addrin = "../../data/COMB/";
            string addrout = "../../results/SPMInput/COAD_CCS/";

            // vaf matrix V and binarized matrix M
            vector<vector<int>> M;
            vector<vector<double>> V;

            readMatrix(addrin+in, V);
            readMatrix(addrin+in, M, true);

            // generate matrix file
            writeMatrix(addrout+out, M, V, K, T);
    }

    return 0;
}