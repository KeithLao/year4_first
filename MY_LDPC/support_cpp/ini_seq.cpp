#include "ini_seq.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
using namespace std;

IniSeq::IniSeq(int a, int b, int c){
    n = a;
    m = b;
    k = c;
}


vector<int> IniSeq::generateRandom(int seed) {
    vector<int> seq(k);
    mt19937 rng(seed);
    uniform_int_distribution<int> dist(0, 1);
    
    for (int i = 0; i < k; i++) {
        seq[i] = dist(rng);  // Generate 0 or 1 sequence
    }
    return seq;
}

vector<int> IniSeq::Encoder(const BaseMatrix& H, const vector<int>& ini){
    using u64 = uint64_t;
    
    if ((int)ini.size() != k) throw std::invalid_argument("ini.size() must equal k = n - m");
    if ((int)H.tanner_graph_i.size() != m) throw std::invalid_argument("H.tanner_graph_i size must equal H.m");

    // number of 64-bit words needed to hold m rows
    int W_rows = (m + 63) / 64;
    // number of 64-bit words needed to hold m columns when solving Hp (for row representation)
    int W_cols = (m + 63) / 64;

    // Build column bitsets: cols[v] is bitset of length m (rows) describing H[:, v]
    vector<vector<u64>> cols(n, vector<u64>(W_rows, 0));
    for (int r = 0; r < m; ++r) {
        for (int v : H.tanner_graph_i[r]) {
            if (v < 0 || v >= n) continue;
            int word = r >> 6;
            int bit = r & 63;
            cols[v][word] |= (u64(1) << bit);
        }
    }

    // 1) Greedy selection of m independent columns (parity columns) using incremental basis insertion.
    // We'll maintain basis vectors keyed by pivot row index. pivot_basis[p] = vector<u64> if pivot exists
    vector<vector<u64>> pivot_basis(m); // empty means no basis with that pivot yet
    vector<int> parity_cols; parity_cols.reserve(m);

    auto highest_set_bit = [&](const vector<u64>& v)->int{
        // return highest row index (0..m-1) of bit set in v, or -1 if zero.
        for (int wi = W_rows - 1; wi >= 0; --wi) {
            u64 w = v[wi];
            if (w) {
                int lead = 63 - __builtin_clzll(w); // highest bit index in this 64-bit word (0..63)
                return wi * 64 + lead;
            }
        }
        return -1;
    };

    auto xor_inplace = [&](vector<u64>& a, const vector<u64>& b){
        for (int i = 0; i < W_rows; ++i) a[i] ^= b[i];
    };

    for (int col = 0; col < n && (int)parity_cols.size() < m; ++col) {
        vector<u64> vec = cols[col]; // column bitset copy
        // reduce with existing pivot_basis
        int piv;
        while ((piv = highest_set_bit(vec)) != -1) {
            if (!pivot_basis[piv].empty()) {
                xor_inplace(vec, pivot_basis[piv]);
            } else {
                // new independent column with pivot at 'piv'
                pivot_basis[piv] = vec; // store basis vector
                parity_cols.push_back(col);
                break;
            }
        }
        // if vec becomes zero -> dependent, skip
    }
    if ((int)parity_cols.size() != m) throw std::runtime_error("Failed to find m independent parity columns (unexpected)");

    // 2) systematic columns are the rest (in original order)
    vector<char> is_parity(n, 0);
    for (int c : parity_cols) is_parity[c] = 1;
    vector<int> sys_cols; sys_cols.reserve(k);
    for (int c = 0; c < n; ++c) if (!is_parity[c]) sys_cols.push_back(c);
    if ((int)sys_cols.size() != k) throw std::runtime_error("systematic column count mismatch");

    // 3) permutation: perm[new_index] = original_col (sys_cols first, then parity_cols)
    vector<int> perm; perm.reserve(n);
    for (int c : sys_cols) perm.push_back(c);
    for (int c : parity_cols) perm.push_back(c);

    // 4) Build Hp as m x m matrix in row-wise bitset form, and Hs implicitly via cols.
    // Hp_rows[r] is bitset length m representing Hp[r, :]
    vector<vector<u64>> Hp_rows(m, vector<u64>(W_cols, 0));
    // For each parity column j (0..m-1), parity_col = perm[k + j]
    for (int j = 0; j < m; ++j) {
        int orig_col = perm[k + j];
        // column bits of this parity column are cols[orig_col] (length m bits)
        // for each row r where bit is set, set bit j in Hp_rows[r]
        for (int wi = 0; wi < W_rows; ++wi) {
            u64 word = cols[orig_col][wi];
            while (word) {
                int tz = __builtin_ctzll(word);
                int r = (wi << 6) + tz;
                if (r < m) {
                    int wj = j >> 6;
                    int bj = j & 63;
                    Hp_rows[r][wj] |= (u64(1) << bj);
                }
                word &= word - 1;
            }
        }
    }

    // 5) Build s_perm from ini (ini is assumed in correct permuted systematic order)
    vector<unsigned char> s_perm(k, 0);
    for (int i = 0; i < k; ++i) s_perm[i] = (unsigned char)(ini[i] & 1);

    // 6) Compute lambda = Hs * s_perm  (Hs is columns perm[0..k-1])
    vector<unsigned char> lambda(m, 0);
    for (int i = 0; i < k; ++i) {
        if (!s_perm[i]) continue;
        int orig_col = perm[i];
        // for each row r where cols[orig_col] has 1, flip lambda[r]
        for (int wi = 0; wi < W_rows; ++wi) {
            u64 word = cols[orig_col][wi];
            while (word) {
                int tz = __builtin_ctzll(word);
                int r = (wi << 6) + tz;
                if (r < m) lambda[r] ^= 1;
                word &= word - 1;
            }
        }
    }

    // 7) Solve Hp * p = lambda  using GF(2) Gaussian elimination on rows (in-place)
    // We'll perform elimination on Hp_rows with rhs lambda.
    // Make local copy so we can modify
    vector<vector<u64>> A = Hp_rows; // m rows, each has W_cols words
    vector<unsigned char> b = lambda;

    // forward elimination
    for (int col = 0, row = 0; col < m && row < m; ++col) {
        // find pivot row with bit at 'col' set, starting from row
        int pivot = -1;
        int widx = col >> 6;
        u64 mask = u64(1) << (col & 63);
        for (int r = row; r < m; ++r) {
            if (A[r][widx] & mask) { pivot = r; break; }
        }
        if (pivot == -1) continue; // no pivot in this column -> should not happen (we chose independent cols)
        if (pivot != row) {
            swap(A[pivot], A[row]);
            swap(b[pivot], b[row]);
        }
        // eliminate below
        for (int r = row + 1; r < m; ++r) {
            if (A[r][widx] & mask) {
                // row r ^= row
                for (int wi = 0; wi < W_cols; ++wi) A[r][wi] ^= A[row][wi];
                b[r] ^= b[row];
            }
        }
        ++row;
    }

    // back substitution to obtain p (length m)
    vector<unsigned char> p(m, 0);
    // find pivot positions per row (leading one)
    vector<int> lead_col(m, -1);
    for (int r = 0; r < m; ++r) {
        // find first set bit in row r
        int found = -1;
        for (int wi = 0; wi < W_cols; ++wi) {
            u64 word = A[r][wi];
            if (word) {
                int tz = __builtin_ctzll(word);
                found = (wi << 6) + tz;
                break;
            }
        }
        lead_col[r] = found; // could be -1 if zero row
    }

    for (int r = m - 1; r >= 0; --r) {
        int lc = lead_col[r];
        if (lc == -1) {
            if (b[r]) throw std::runtime_error("No solution in Hp * p = lambda (inconsistent)"); // shouldn't happen
            continue;
        }
        // compute sum of A[r][c] * p[c] for c > lc
        unsigned char acc = b[r];
        // iterate over bits in A[r]; faster to iterate words
        for (int wi = (lc >> 6); wi < W_cols; ++wi) {
            u64 word = A[r][wi];
            // mask out bits <= lc
            if (wi == (lc >> 6)) {
                int bitpos = lc & 63;
                u64 masklow = ((u64(1) << (bitpos + 1)) - 1);
                word &= ~masklow;
            }
            while (word) {
                int tz = __builtin_ctzll(word);
                int c = (wi << 6) + tz;
                if (c > lc && c < m) acc ^= p[c];
                word &= word - 1;
            }
        }
        p[lc] = acc;
    }

    // 8) assemble x_perm = [s_perm | p]
    vector<int> x_perm(n, 0);
    for (int i = 0; i < k; ++i) x_perm[i] = s_perm[i];
    for (int i = 0; i < m; ++i) x_perm[k + i] = p[i];

    // 9) inverse permute to original ordering
    vector<int> codeword(n, 0);
    for (int newidx = 0; newidx < n; ++newidx) {
        int orig_col = perm[newidx];
        codeword[orig_col] = x_perm[newidx];
    }

    // 10) verify parity-check (optional but kept in debug builds). We keep but it's cheap.
    for (int r = 0; r < m; ++r) {
        unsigned char acc = 0;
        for (int v : H.tanner_graph_i[r]) acc ^= (unsigned char)(codeword[v] & 1);
        if (acc != 0) throw std::runtime_error("Parity-check failed after encoding");
    }

    // 11) map (0,1) to (1,-1)
    vector<int> messages(n);
    for (int i = 0; i < n; ++i) messages[i] = 1 - (2 * codeword[i]);
    return messages;
}

vector<double> IniSeq::addAWGN(const vector<int>& input, double snr_db, int seed, double bitrate) {
    vector<double> output(n);
    const double signal_power = 1.0;                          // ±1 sequence => Es = 1
    double snr_linear = std::pow(10.0, snr_db / 10.0);        // dB -> linear

    // real AWGN for BPSK: sigma^2 = Es / (2 * (Es/N0)) = signal_power / (2 * snr_linear)
    double noise_variance = (snr_linear > 0.0)
                            ? (signal_power / (2.0 * snr_linear * bitrate))
                            : std::numeric_limits<double>::infinity();
    double noise_std = std::sqrt(noise_variance);

    std::mt19937 rng(seed);
    std::normal_distribution<double> noise_dist(0.0, noise_std);
    for (size_t i = 0; i < n; ++i) output[i] = static_cast<double>(input[i]) + noise_dist(rng);

    return output;
}

vector<double> IniSeq::generate_noise(double snr) {
    vector<double> output(n);
    // RV generator
    static thread_local std::mt19937_64 rng((std::random_device())());

    double snr_linear = std::pow(10.0, snr / 10.0); // snr (linear) = Eb/N0
    double Eb = 1.0;                                // (1,-1)seq. -> Eb = 1.0
    double N0 = Eb / snr_linear;                    // N0 = Eb / (Eb/N0)
    double eta = 1.0;                               // eta = 1.0 for peg508x1004
    double sigma2 = (eta * eta) * N0 / 2.0;         // sigma^2 = η^2 * N0 / 2
    double sigma = (sigma2 > 0.0) ? std::sqrt(sigma2) : 0.0;

    std::normal_distribution<double> stdn(0.0, sigma);
    for(int i = 0; i < n; i++){
        double noise = stdn(rng);
        output[i] = noise;
        //cout << "Noise value: "<< i <<" th " << output[i] << endl;
    }
    //cout << "SNR: "<< snr <<" sigma: " << sigma << endl;
    return output;
}

vector<int> IniSeq::hard_decision(const vector<double>& received){
    vector<int> hard_received(n);
    for (int i = 0; i < n; i++) {
        hard_received[i] = (received[i] >= 0) ? 1 : -1;
    }
    return hard_received;
}