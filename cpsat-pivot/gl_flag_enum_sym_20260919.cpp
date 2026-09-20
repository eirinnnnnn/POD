// gl_flag_enum_sym -- exact enumeration of pivot profiles over the GL(m,2) orbit of C' = C_b P_field, m <= 8, k <= 128.
//
//   profile(A) = profile(A L)          L = lower-unitriangular (transvections z_i<-z_i+z_j, i>j are profile-neutral)
//   profile(A) = profile(gamma A)      gamma in Gamma = stabiliser of the code (Singer cycle x->alpha x, Frobenius)
//   Gamma contains the Singer cycle, which is TRANSITIVE on nonzero vectors, so the last column of A may be fixed to 1:
//   flags to visit = [m]_2!/(2^m-1)   (m=7: 615,195   m=8: 78,129,765)      [fixlast=1]
//
// usage: gl_flag_enum_sym cp.txt [fixlast=0|1] [maxprefix=0]      (maxprefix>0: benchmark on the first prefixes only)
// cp.txt: "m k" then 2^m lines "hi lo" (hex 64-bit) = column x of G_b[:,pi_field] (k bits).
// stdout: #flags, then per distinct profile:  count A_0..A_{m-1} profile_hex(64 hex chars, bit j = index j)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <omp.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;
static int M, N, K;
static u128 CP[256];
struct Key { u64 w[4]; bool operator==(const Key& o) const { return w[0]==o.w[0]&&w[1]==o.w[1]&&w[2]==o.w[2]&&w[3]==o.w[3]; } };
struct KHash { size_t operator()(const Key& k) const { u64 h = 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 4; i++) h ^= k.w[i] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); return (size_t)h; } };
struct Rec { u64 count; int cols[8]; };
typedef std::unordered_map<Key, Rec, KHash> Map;
static inline int hb128(u128 x) { u64 hi = (u64)(x >> 64); if (hi) return 127 - __builtin_clzll(hi); return 63 - __builtin_clzll((u64)x); }
static Key profile_of(const int* A) {
    int Q[256]; Q[0] = 0;
    for (int a = 1; a < N; a++) Q[a] = Q[a & (a - 1)] ^ A[__builtin_ctz(a)];
    u128 w[256];
    for (int a = 0; a < N; a++) w[a] = CP[Q[a]];
    for (int i = 0; i < M; i++) { int b = 1 << i; for (int x = 0; x < N; x++) if (!(x & b)) w[x] ^= w[x | b]; }
    u128 basis[128]; memset(basis, 0, sizeof(basis));
    Key key; memset(&key, 0, sizeof(key));
    for (int j = 0; j < N; j++) {
        u128 v = w[j];
        while (v) { int h = hb128(v);
            if (!basis[h]) { basis[h] = v; key.w[j >> 6] |= 1ULL << (j & 63); break; }
            v ^= basis[h]; }
    }
    return key;
}
struct State { int A[8]; int pivmask; };
static int FIXLAST = 0;
static void dfs(State& s, int depth, Map& mp, u64& leaves) {
    if (depth == M) {
        Key p = profile_of(s.A);
        auto it = mp.find(p);
        if (it == mp.end()) { Rec r; r.count = 1; for (int i = 0; i < M; i++) r.cols[i] = s.A[i]; mp.emplace(p, r); }
        else it->second.count++;
        leaves++; return;
    }
    int c = M - 1 - depth;
    for (int v = 1; v < N; v++) {
        if (v & s.pivmask) continue;
        State t = s; t.A[c] = v; t.pivmask |= (1 << (31 - __builtin_clz(v)));
        dfs(t, depth + 1, mp, leaves);
    }
}
static void prefixes(State& s, int depth, int stop, std::vector<State>& out) {
    if (depth == stop) { out.push_back(s); return; }
    int c = M - 1 - depth;
    for (int v = 1; v < N; v++) {
        if (depth == 0 && FIXLAST && v != 1) continue;   // Singer-cycle symmetry: last column := 1
        if (v & s.pivmask) continue;
        State t = s; t.A[c] = v; t.pivmask |= (1 << (31 - __builtin_clz(v)));
        prefixes(t, depth + 1, stop, out);
    }
}
int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s cp.txt [fixlast] [maxprefix]\n", argv[0]); return 1; }
    FIXLAST = argc > 2 ? atoi(argv[2]) : 0; long maxp = argc > 3 ? atol(argv[3]) : 0;
    FILE* f = fopen(argv[1], "r"); if (!f || fscanf(f, "%d %d", &M, &K) != 2) return 2;
    if (K > 128 || M > 8) { fprintf(stderr, "need k<=128, m<=8\n"); return 4; }
    N = 1 << M;
    for (int x = 0; x < N; x++) { unsigned long long hi, lo; if (fscanf(f, "%llx %llx", &hi, &lo) != 2) return 3; CP[x] = ((u128)hi << 64) | lo; }
    fclose(f);
    State s0; memset(&s0, 0, sizeof(s0));
    std::vector<State> pre; prefixes(s0, 0, FIXLAST ? 3 : 2, pre);
    long np = (long)pre.size(); if (maxp > 0 && maxp < np) np = maxp;
    int T = omp_get_max_threads(); std::vector<Map> maps(T); std::vector<u64> leaves(T, 0);
    fprintf(stderr, "m=%d k=%d fixlast=%d prefixes=%zu (running %ld) threads=%d\n", M, K, FIXLAST, pre.size(), np, T);
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < np; i++) { int t = omp_get_thread_num(); State s = pre[i]; dfs(s, FIXLAST ? 3 : 2, maps[t], leaves[t]); }
    double dt = omp_get_wtime() - t0;
    Map all; u64 total = 0;
    for (int t = 0; t < T; t++) { total += leaves[t];
        for (auto& kv : maps[t]) { auto it = all.find(kv.first); if (it == all.end()) all.emplace(kv.first, kv.second); else it->second.count += kv.second.count; } }
    fprintf(stderr, "flags=%llu  time=%.2fs  -> %.2f us/flag/thread-equiv %.2f us/flag wall\n", total, dt, dt * T * 1e6 / (double)total, dt * 1e6 / (double)total);
    printf("%llu\n", total);
    for (auto& kv : all) { printf("%llu", kv.second.count); for (int i = 0; i < M; i++) printf(" %d", kv.second.cols[i]);
        printf(" %016llx%016llx%016llx%016llx\n", kv.first.w[3], kv.first.w[2], kv.first.w[1], kv.first.w[0]); }
    return 0;
}
