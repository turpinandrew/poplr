#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <set>
#include <boost/container_hash/hash.hpp>
#include <chrono>

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
// Always have 1..n as the first (not so random)
class RandomPermsIterator: public PermutationBase {
public:
    RandomPermsIterator(int n, int count) : PermutationBase(n, count) {
        std::srand(std::chrono::system_clock::now().time_since_epoch().count());
        ks.insert(0);
        while (ks.size() < count)
            ks.insert(rand() % fact[n]);
        kIt = ks.begin();
    }
    using PermutationBase::PermutationBase;

    bool hasNext() const override {
        return kIt != ks.end();
    }

    std::vector<int> next() override {
        if (kIt == ks.end())
            return std::vector<int>();

        int k = *kIt++;
        return getPermutation(n, k);
    }

    // Get the k'th permutation of 1..n (k >= 0)
    // Adaption of  Stackoverflow answer: Aug 21, 2014 at 21:33 Ismael EL ATIFI
    std::vector<int> getPermutation(int n, long k) {
        std::set<int> nums;
        for(int i = 1 ; i <= n ; i++)
            nums.insert(i);
        
        std::vector<int> result;

        while (nums.size() > 1) {
            long sizeGroup = fact[nums.size() - 1];
            
            long q = k / sizeGroup;
            long r = k - q * sizeGroup;

            int val = *std::next(nums.begin(), q);
            
            result.push_back(val);  
            nums.erase(val);

            k = r;
        }

        result.push_back(*nums.begin());  // the final digit remaining in nums
        
        return result;
    }

private:
    long fact[15] = {1,1,2,6,24,120,720,5040,40320,362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};

    std::unordered_set<int> ks;
    std::unordered_set<int>::iterator kIt = ks.begin();
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