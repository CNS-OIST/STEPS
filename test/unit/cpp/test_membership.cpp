#include <vector>

#include <algorithm>
#include <vector>
#include <set>
#include <list>

#include "steps/util/collections.hpp"

#include "gtest/gtest.h"

using namespace steps::util;

// global arrays for test data

float X[]={5,3,1,2,4,6,5,7,2,3};
float E[]={10,9,8,7,6,5,4,3,2,1,0};

#define ARRSZ(a) (sizeof(a)/sizeof((a)[0]))
constexpr int n_X=ARRSZ(X);
constexpr int n_E=ARRSZ(E);

bool   M[n_E]={false,false,false,true,true,true,true,true,true,true,false};


TEST(Membership,SimpleScalar) {
    auto v=map_membership(E,X);
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));
}

TEST(Membership,MixedContainers) {
    std::vector<float> x_vec(std::begin(X),std::end(X));
    std::list<float>   e_list(std::begin(E),std::end(E));

    auto v=map_membership(e_list,x_vec);
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));
}

template <typename V>
struct count_copies {
    count_copies() =default;

    count_copies(V v_): v(v_) {}
    count_copies(const count_copies &x): v(x.v) { ++copies; }
    count_copies &operator=(const count_copies &x) {
        if (this!=&x) v=x.v, ++copies;
        return *this;
    }

    operator V() const { return v; } 
    V v;

    static int copies;
};

template <typename V> int count_copies<V>::copies=0;

TEST(Membership,ConfirmCopies) {
    count_copies<float> cc_x[n_X];
    count_copies<float> cc_e[n_E];

    std::copy(std::begin(X),std::end(X),std::begin(cc_x));
    std::copy(std::begin(E),std::end(E),std::begin(cc_e));

    count_copies<float>::copies=0;

    // With the hash_references tag, no copies should be made of elements of x.
    hash_references_tag hash_references;

    auto v=map_membership(hash_references,cc_e,cc_x,std::hash<double>());
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));
    ASSERT_EQ(count_copies<float>::copies,0);

    // Without the hash_references tag, the internal hash table will have to
    // copy at least the unique items in x

    size_t n_X_unique=std::set<float>(std::begin(X),std::end(X)).size();

    v=map_membership(cc_e,cc_x,std::hash<double>());
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));

    ASSERT_GE(count_copies<float>::copies,n_X_unique);
}

TEST(Membership,ImplicitConversion) {
    double dbl_x[n_X];
    std::copy(std::begin(X),std::end(X),std::begin(dbl_x));

    auto v=map_membership(E,dbl_x);
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));
}

TEST(Membership,CustomHash) {
    int silly_hash_count=0;

    auto v=map_membership(E,X,[&silly_hash_count](float x) -> size_t { return ++silly_hash_count,0; });
    ASSERT_TRUE(std::equal(v.begin(),v.end(),std::begin(M)));

    ASSERT_GT(silly_hash_count,0);
}
