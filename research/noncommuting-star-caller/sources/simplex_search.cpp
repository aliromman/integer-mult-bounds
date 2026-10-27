// Apache-2.0. Ali Romman, with substantial OpenAI Codex assistance.
// The simplex and projector parametrization are due to Cinna Davis.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using Mat = uint16_t;
constexpr Mat ID = 0x8421;
constexpr int N = 4;
constexpr int L = 5;
using Poly = std::array<Mat, L + 1>;

Mat mul(Mat a, Mat b) {
    Mat out = 0;
    for (int i = 0; i < N; ++i) {
        int row = 0;
        for (int j = 0; j < N; ++j)
            if ((a >> (4*i+j)) & 1) row ^= (b >> (4*j)) & 15;
        out |= row << (4*i);
    }
    return out;
}

Mat outer(int u, int v) {
    Mat out = 0;
    for (int i = 0; i < N; ++i)
        if ((u >> i) & 1) out |= v << (4*i);
    return out;
}

int rank(Mat a) {
    int pivots[4] = {}, count = 0;
    for (int i = 0; i < N; ++i) {
        int row = (a >> (4*i)) & 15;
        for (int k = 3; k >= 0; --k) {
            if (!((row >> k) & 1)) continue;
            if (pivots[k]) row ^= pivots[k];
            else { pivots[k] = row; ++count; break; }
        }
    }
    return count;
}

int parity(int x) { return __builtin_parity(static_cast<unsigned>(x)); }

std::array<Mat, L> r;
std::array<std::vector<Mat>, L> ps;
std::array<Mat, L> word;
uint64_t total = 0, all_query = 0, factored = 0;
std::map<int, uint64_t> rank_histogram;
std::map<std::array<int, L+1>, uint64_t> degree_histogram;
int full_mask = 31;
int inverse_mask = 0;

void visit(int depth, const Poly& a, const Poly& b) {
    if (depth == L) {
        ++total;
        if (a != b) return;
        ++all_query;
        std::array<int, L+1> ranks;
        int sum = 0;
        for (int d = 0; d <= L; ++d) { ranks[d] = rank(a[d]); sum += ranks[d]; }
        ++rank_histogram[sum];
        if (sum != N) return;
        ++factored;
        ++degree_histogram[ranks];
        if (factored <= 10) {
            std::cout << "WORD";
            for (Mat p : word) std::cout << ' ' << p;
            std::cout << " COEFF";
            for (Mat c : a) std::cout << ' ' << c;
            std::cout << '\n';
        }
        return;
    }
    for (Mat p : ps[depth]) {
        word[depth] = p;
        Poly aa{}, bb{};
        for (int d = 0; d <= depth; ++d) {
            Mat pa = mul(p, a[d]), pb = mul(p, b[d]);
            bool inverse = (inverse_mask >> depth) & 1;
            aa[d] ^= inverse ? pa : a[d] ^ pa;
            aa[d+1] ^= inverse ? a[d] ^ pa : pa;
            bb[d] ^= inverse ? pb : b[d] ^ pb;
            bb[d+1] ^= inverse ? b[d] ^ pb : pb;
            bb[d + (((full_mask >> depth) & 1) && !inverse)] ^= mul(r[depth], a[d]);
        }
        visit(depth+1, aa, bb);
    }
}

int main(int argc, char** argv) {
    if (argc >= 2) full_mask = std::atoi(argv[1]);
    if (argc >= 3) inverse_mask = std::atoi(argv[2]);
    if ((inverse_mask & full_mask) != inverse_mask) return 2;
    for (int i = 0; i < L; ++i) {
        int u = i < N ? 1 << i : 15;
        int v = i < N ? 15 ^ u : 15;
        r[i] = outer(u, v);
        std::set<Mat> unique;
        for (int a = 1; a < 16; ++a) {
            if (!parity(v & a)) continue;
            for (int b = 1; b < 16; ++b) {
                if (!parity(u & b) || parity(a & b)) continue;
                unique.insert(outer(u, b) ^ outer(a, v));
            }
        }
        ps[i].assign(unique.begin(), unique.end());
        if (!((full_mask >> i) & 1)) ps[i] = {0};
        std::cout << "INVENTORY " << i << ' ' << r[i] << ' ' << ps[i].size() << '\n';
        for (Mat p : ps[i]) {
            if (!p) continue;
            if (rank(p) != 2 || mul(p,p) != p || mul(p,r[i]) != r[i] || mul(r[i],p) != r[i])
                return 2;
        }
    }
    Poly a{}, b{};
    a[0] = ID;
    visit(0, a, b);
    std::cout << "TOTAL " << total << " ALL_QUERY " << all_query << " FACTORED " << factored << '\n';
    for (const auto& [k, v] : rank_histogram) std::cout << "RANK_SUM " << k << ' ' << v << '\n';
    for (const auto& [k, v] : degree_histogram) {
        std::cout << "DEGREES";
        for (int x : k) std::cout << ' ' << x;
        std::cout << " COUNT " << v << '\n';
    }
}
