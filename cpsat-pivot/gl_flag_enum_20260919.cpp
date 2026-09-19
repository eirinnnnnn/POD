// gl_flag_enum -- EXACT enumeration of the pivot profiles reachable in the GL(m,2) orbit of C' = C_b P_field.
//
// Fact (verified numerically, m=6 and m=7): appending a transvection T=(i,j), z_i <- z_i + z_j with i>j, on the
// RIGHT of A never changes the column pivot profile.  Those T generate the lower-unitriangular group L, so
//        profile(A) = profile(A L)   for all L,
// and the profile is a function of the coset A*L, i.e. of the FLAG  span(A_{m-1}) < span(A_{m-2},A_{m-1}) < ...
// There are prod_{c=1..m}(2^c-1) = [m]_2! flags (78,129,765 for m=7) instead of |GL(m,2)| = 1.6e14.
//
// Per candidate A (columns A_i = A e_i as ints, bit 0 = LSB):
//   Q_A(a) = A a ;  newcol[a] = CP[Q_A(a)]   (CP[x] = column x of G_b[:, pi_field], k bits)
//   W[:,j] = XOR_{a superset of j} newcol[a]           (superset-sum = F^(x m) on columns)
//   profile = columns of W that are independent of all earlier columns
//
// Input file:  "m k" then 2^m lines "hi lo" (64-bit hex each) = CP[x].
// Output:      first line = number of flags visited (must equal [m]_2!), then one line per DISTINCT profile:
//              count  A_0 ... A_{m-1}  profile_hex
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <omp.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;

static int M, N, K;
static u128 CP[128];

struct Rec { u64 count; int cols[8]; };
struct Hash { size_t operator()(const u128& x) const {
    u64 a = (u64)x, b = (u64)(x >> 64);
    return (size_t)((a ^ (b * 0x9E3779B97F4A7C15ULL)) * 0xBF58476D1CE4E5B9ULL); } };
typedef std::unordered_map<u128, Rec, Hash> Map;

static inline int hb128(u128 x) {
    u64 hi = (u64)(x >> 64);
    if (hi) return 127 - __builtin_clzll(hi);
    return 63 - __builtin_clzll((u64)x);
}

static u128 profile_of(const int* A) {
    int Q[128]; Q[0] = 0;
    for (int a = 1; a < N; a++) Q[a] = Q[a & (a - 1)] ^ A[__builtin_ctz(a)];
    u128 w[128];
    for (int a = 0; a < N; a++) w[a] = CP[Q[a]];
    for (int i = 0; i < M; i++) {
        int b = 1 << i;
        for (int x = 0; x < N; x++) if (!(x & b)) w[x] ^= w[x | b];
    }
    u128 basis[128]; memset(basis, 0, sizeof(basis));
    u128 prof = 0;
    for (int j = 0; j < N; j++) {
        u128 v = w[j];
        while (v) {
            int h = hb128(v);
            if (!basis[h]) { basis[h] = v; prof |= ((u128)1 << j); break; }
            v ^= basis[h];
        }
    }
    return prof;
}

struct State { int A[8]; int pivmask; };

static void dfs(State& s, int depth, Map& mp, u64& leaves) {
    if (depth == M) {
        u128 p = profile_of(s.A);
        auto it = mp.find(p);
        if (it == mp.end()) { Rec r; r.count = 1; for (int i = 0; i < M; i++) r.cols[i] = s.A[i]; mp.emplace(p, r); }
        else it->second.count++;
        leaves++;
        return;
    }
    int c = M - 1 - depth;                       // columns are fixed from the LAST one backwards
    for (int v = 1; v < N; v++) {
        if (v & s.pivmask) continue;             // canonical coset representative: zero on all pivot bits
        int h = 31 - __builtin_clz(v);
        State t = s; t.A[c] = v; t.pivmask |= (1 << h);
        dfs(t, depth + 1, mp, leaves);
    }
}

static void prefixes(State& s, int depth, int stop, std::vector<State>& out) {
    if (depth == stop) { out.push_back(s); return; }
    int c = M - 1 - depth;
    for (int v = 1; v < N; v++) {
        if (v & s.pivmask) continue;
        int h = 31 - __builtin_clz(v);
        State t = s; t.A[c] = v; t.pivmask |= (1 << h);
        prefixes(t, depth + 1, stop, out);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s cp.txt\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "r");
    if (!f || fscanf(f, "%d %d", &M, &K) != 2) return 2;
    N = 1 << M;
    for (int x = 0; x < N; x++) {
        unsigned long long hi, lo;
        if (fscanf(f, "%llx %llx", &hi, &lo) != 2) return 3;
        CP[x] = ((u128)hi << 64) | lo;
    }
    fclose(f);
    State s0; memset(&s0, 0, sizeof(s0));
    std::vector<State> pre; prefixes(s0, 0, 2, pre);
    int T = omp_get_max_threads();
    std::vector<Map> maps(T); std::vector<u64> leaves(T, 0);
    fprintf(stderr, "m=%d k=%d  prefixes=%zu  threads=%d\n", M, K, pre.size(), T);
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < (long)pre.size(); i++) {
        int t = omp_get_thread_num();
        State s = pre[i];
        dfs(s, 2, maps[t], leaves[t]);
    }
    Map all; u64 total = 0;
    for (int t = 0; t < T; t++) {
        total += leaves[t];
        for (auto& kv : maps[t]) {
            auto it = all.find(kv.first);
            if (it == all.end()) all.emplace(kv.first, kv.second);
            else it->second.count += kv.second.count;
        }
    }
    printf("%llu\n", total);
    for (auto& kv : all) {
        printf("%llu", kv.second.count);
        for (int i = 0; i < M; i++) printf(" %d", kv.second.cols[i]);
        printf(" %016llx%016llx\n", (u64)(kv.first >> 64), (u64)kv.first);
    }
    return 0;
}
