#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>

// Base class for permutation iterators
class PermutationBase {
public:
    PermutationBase() { this->n = 0 ; this->count = 0; }
    PermutationBase(int n) { this->n = n ; this->count = n; }
    PermutationBase(int n, int count) { this->n = n ; this->count = count; }
    virtual bool hasNext() const = 0 ;
    virtual std::vector<int> next() = 0;

protected:
    int n, count;
};

// Generate count random permutations of 1..n
class RandomPermsIterator: public PermutationBase {
public:
    RandomPermsIterator(int n, int count) : nums(n), PermutationBase(n, count) {
        for (int i = 0; i < n; i++)
            nums[i] = i + 1;

        perms.insert(nums);
    }
    using PermutationBase::PermutationBase;

    bool hasNext() const override {
        return count > 0;
    }

    std::vector<int> next() override {
        if (count == 0) 
            return std::vector<int>();

        bool inserted;
        std::vector<int> b = nums;
        do {
            std::srand(std::time(0));
            std::shuffle(nums.begin(), nums.end(), std::default_random_engine(std::rand()));
        } while (perms.count(nums));
        perms.insert(nums);

        count--;

        return b;
    }
private:
    std::vector<int> nums;
    std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> perms;
};

// Heap's algorithm for generating all permutations of 1..n
// adapted from Wikipedia as an iterator
class AllPermsIterator: public PermutationBase {
public:
    AllPermsIterator(int n) : A(n), c(n), PermutationBase(n) {
        for (int j = 0; j < n; j++)
            A[j] = j + 1;

        for (int j = 0; j < n; j++)
            c[j] = 0;

        i = 1; 
    }
    using PermutationBase::PermutationBase;

    std::vector<int> next() override {
        if (i == n)
            return std::vector<int>();

        std::vector<int> output = A;

        while (c[i] >= i && i < n) {
            c[i] = 0;
            i++;
        }

        if (i < n) {
            if (i % 2 == 0)
                std::swap(A[0], A[i]);
            else
                std::swap(A[c[i]], A[i]);
            c[i]++;
            i = 1;
        }

        return output;
    }

    bool hasNext() const override { return i < n; }
    
private:
    std::vector<int> A;
    std::vector<int> c;
    int i;
    std::vector<int> output;
};

// Iterator for all permutations of 1..n
// Chooses correct iterator based on count relative to n
class PermutationIterator {
public:
    PermutationIterator(int n, int count) {
        if (count >= std::tgamma(n + 1)) {
            api = AllPermsIterator(n);
            it = &api;
        } else {
            rpi = RandomPermsIterator(n, count);
            it = &rpi;
        }
    }

    bool hasNext() const {return it->hasNext();}

    std::vector<int> next() const {return it->next();}

private:
    PermutationBase *it;
    AllPermsIterator api;
    RandomPermsIterator rpi;
};