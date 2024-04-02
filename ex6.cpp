#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <iterator>

#include <cmath>
#include <vector>
#include <version>

#ifdef NOFORMAT
#define HAVEFORMAT 0
#else
#ifdef __cpp_lib_format
#define HAVEFORMAT 1
#endif
#endif

using std::pair, std::vector;

struct dvec2 : public pair<double, double> {
    // inherit constructors
    using pair<double, double>::pair;

    // returns the euclidean 2-norm
    double norm() const { return std::sqrt(dot(*this, *this)); }

    // returns the dot product of a and b.
    friend double dot(const dvec2& a, const dvec2& b)
    {
        return a.first * b.first + a.second * b.second;
    }

    // returns the angle theta between a and b.
    friend double theta(const dvec2& a, const dvec2& b)
    {
        return std::acos(dot(a, b) / (a.norm() * b.norm()));
    }

    friend std::ostream& operator<<(std::ostream& os, const dvec2& self)
    {
        return os << "[" << self.first << ", " << self.second << "]";
    }

    friend bool operator==(const dvec2&, const dvec2&) = default;
};

#if HAVEFORMAT
#include <format>
template<>
struct std::formatter<dvec2> : std::formatter<std::string> {
    template<typename FmtCtx>
    auto format(dvec2 v, FmtCtx& ctx) const
    {
        return std::formatter<std::string>::format(
            std::format("[{:g}, {:g}]", v.first, v.second), ctx);
    }
};
#endif

vector<dvec2> ingest_dvecs(std::istream& input)
{
    std::istream_iterator<double> instream{input}, end;
    vector<dvec2> output;

    while (instream != end) {
        auto x = *instream++;
        if (instream == end)
            throw std::runtime_error("mismatched vector elements");
        auto y = *instream++;
        output.push_back({x, y});
    }
    return output;
}

// returns all unique (i.e. [a,b]==[b,a]) pairs of vectors in `vecs`, excluding
// pairs of the same vector.
vector<pair<dvec2, dvec2>> pairwise_elts(const vector<dvec2>& vecs)
{
    vector<pair<dvec2, dvec2>> pairs;

    auto it = vecs.cbegin();
    auto end = vecs.cend();
    for (; it != end; ++it) {
        for (auto it2 = it + 1; it2 != end; ++it2) {
            pairs.push_back({*it, *it2});
        }
    }
    return pairs;
}

// return the pairs of dvec2s ordered by theta in ascending order
vector<pair<dvec2, dvec2>> theta_sort(const vector<dvec2>& vecs)
{
    auto pairs = pairwise_elts(vecs);
    std::sort(pairs.begin(), pairs.end(), [](auto x, auto y) {
        return theta(x.first, x.second) < theta(y.first, y.second);
    });
    return pairs;
}

#ifndef TESTING
const std::string DEFAULT_FNAME = "test.txt";

int main(int argc, char** argv)
{
    std::string fname = DEFAULT_FNAME;
    if (argc == 2) {
        fname = argv[1];
    }

    std::ifstream ifile{fname};
    if (!ifile.is_open()) {
        throw std::runtime_error("no input file");
    }

    auto vecs = ingest_dvecs(ifile);
    auto vecpairs = theta_sort(vecs);

    for (auto [x, y] : vecpairs) {
#if HAVEFORMAT
        std::format_to(std::ostreambuf_iterator{std::cout},
                       "𝜃({}, {}) = {:f}\n", x, y, theta(x, y));
#else
        std::cout << "𝜃(" << x << ", " << y << ") = " << std::fixed
                  << theta(x, y) << std::endl
                  << std::defaultfloat;
#endif
    }
    return 0;
}
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

TEST_CASE("ingest_dvecs")
{
    std::istringstream input{"1 1\n1 2\n1 3\n1 4\n1 5"};
    vector<dvec2> expect = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};

    auto result = ingest_dvecs(input);

    CHECK(result == expect);
}

TEST_CASE("theta")
{
    // expected results calculated in Mathematica
    vector<dvec2> input = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
    vector<double> expect = {0.321751, 0.463648, 0.54042,  0.588003,
                             0.141897, 0.218669, 0.266252, 0.0767719,
                             0.124355, 0.0475831};

    auto pairs = pairwise_elts(input);

    for (int i{0}; auto p : pairs) {
        CHECK(theta(p.first, p.second) == doctest::Approx(expect[i++]));
    }
}

TEST_CASE("theta_sort")
{
    vector<dvec2> input = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
    auto result = theta_sort(input);

    // ensure ascending order.
    for (double last{0.0}; auto p : result) {
        auto curr = theta(p.first, p.second);
        CHECK(last <= curr);
        last = curr;
    }
}

TEST_CASE("pairwise_elts")
{
    // just ensuring it returns the correct number of elts and that none are
    // identical.
    vector<dvec2> input = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
    auto result = pairwise_elts(input);
    CHECK(result.size() == 10); // Binomial(5,2) == 10
    for (auto [x, y] : result) {
        CHECK(x != y);
    }
}

#endif
