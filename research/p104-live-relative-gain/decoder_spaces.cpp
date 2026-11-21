#include <array>
#include <bit>
#include <cstdint>
#include <iostream>

using Label = std::uint32_t;
using Inputs = std::array<Label, 111>;

unsigned tensor_rank(const Inputs& roots, const Inputs& labels) {
    std::array<Label, 24> basis{};
    unsigned rank = 0;
    for (std::size_t i = 0; i < roots.size(); ++i) {
        const Label root = roots[i];
        const Label label = labels[i];
        Label row = 0;
        for (unsigned j = 0; j < 3; ++j) {
            if ((label >> j) & 1u) row |= root << (8 * j);
        }
        while (row) {
            const unsigned pivot = std::bit_width(row) - 1;
            if (basis[pivot]) {
                row ^= basis[pivot];
            } else {
                basis[pivot] = row;
                ++rank;
                break;
            }
        }
        if (rank == 21) break;
    }
    return rank;
}

int main() {
    Inputs signatures{};
    std::array<Inputs, 256> roots{};
    std::array<Label, 512> affine{};
    for (auto& x : signatures) if (!(std::cin >> x) || x >= 512) return 1;
    for (auto& lattice : roots) {
        for (auto& x : lattice) if (!(std::cin >> x) || !x || x >= 256) return 1;
    }
    for (Label bit = 1; bit < 512; bit *= 2) {
        Label value;
        if (!(std::cin >> value) || value >= 512) return 1;
        for (Label t = bit; t < 2 * bit; ++t) affine[t] = affine[t ^ bit] ^ value;
    }
    Label extra;
    if (std::cin >> extra) return 1;
    std::array<unsigned, 9> first_success{};
    unsigned total = 0;
    for (Label a = 1; a < 512; ++a) {
        const unsigned p = std::bit_width(a) - 1;
        for (Label b = 1u << (p + 1); b < 512; ++b) {
            if ((b >> p) & 1u) continue;
            const unsigned q = std::bit_width(b) - 1;
            for (Label c = 1u << (q + 1); c < 512; ++c) {
                if (((c >> p) | (c >> q)) & 1u) continue;
                ++total;
                const std::array<Label, 8> members{0, a, b, a ^ b, c, a ^ c, b ^ c, a ^ b ^ c};
                Inputs labels{};
                for (std::size_t i = 0; i < labels.size(); ++i) {
                    for (unsigned j = 0; j < 3; ++j) {
                        labels[i] |= (std::popcount(signatures[i] & members[1u << j]) & 1u) << j;
                    }
                }
                unsigned tries = 0;
                bool success = false;
                for (Label t : members) {
                    if (t && affine[t] <= 256) continue;
                    ++tries;
                    if (tensor_rank(roots[affine[t] & 255u], labels) >= 21) {
                        success = true;
                        ++first_success[tries];
                        break;
                    }
                }
                if (!success) {
                    std::cerr << "missing triple " << a << ' ' << b << ' ' << c << '\n';
                    return 2;
                }
            }
        }
    }
    const std::array<unsigned, 9> expected{0, 107950, 605441, 71850, 2711, 83, 0, 0, 0};
    if (total != 788035 || first_success != expected) return 3;
    std::cout << "three-spaces " << total << "\nminimum certified rank 21\nfirst successes";
    for (unsigned i = 1; i <= 5; ++i) std::cout << ' ' << first_success[i];
    std::cout << '\n';
}
