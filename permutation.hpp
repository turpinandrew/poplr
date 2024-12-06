#ifndef PERMUTATION_HPP
#define PERMUTATION_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <boost/container_hash/hash.hpp>
#include <chrono>

using namespace std;

// Base class for permutation iterators
// These are infinite iterators, so have no end. 
// If they run out of permutations, then they just cycle back to the start.
class PermutationBase {
public:
    /*
        If boundaries are specified as {b1, b2, .., bm} then the permutation of 1..n
        is partitioned into m + 1 blocks: [1, b1], (b1, b2], .., (bm, n]
        
        `count` is then restricted to be no more than max((b[i] - b[i-1])!)

        @param boundaries Right boundaries of NaN blocks (start-1 indexed)
    */
    PermutationBase() { this->count = 0 ; this->n = 0; }
    PermutationBase(int n) {this->count = fact()[n] ; this->n = n; }
    PermutationBase(int n, long count) {this->count = count ; this->n = n; }

    virtual vector<int> next() = 0;

        // based on dyp Dec 31, 2020 at 11:32
        // https://stackoverflow.com/questions/65519616/populating-a-constexpr-array-using-c17
    static constexpr std::array<long, 15> fact() {
        std::array<long, 15> a{1};
        for (auto i = 1; i < a.size(); ++i) {
            a[i] = i * a[i - 1];
        }
        return a;
    }

    virtual ~PermutationBase() = default;

protected:
    int n, count;
};

// Generate count random permutations of 1..n
// Always have 1..n as the first (not so random)
// Does not allow for infinite iteration; returns {} after count iterations.
class RandomPermsIterator: public PermutationBase {
public:
    RandomPermsIterator(int n, long count) : PermutationBase(n, count) {
        assert(count <= fact()[n]);
        srand(chrono::system_clock::now().time_since_epoch().count());
        uniqueKs.insert(0);
        while (uniqueKs.size() < count)
            uniqueKs.insert(rand() % fact()[n]);
        kIt = uniqueKs.begin();
    }
    using PermutationBase::PermutationBase;

    vector<int> next() override {
        int k = *kIt++;
        return getPermutation(n, k);
    }

    // Get the k'th permutation of 1..n (k >= 0)
    // Adaption of  Stackoverflow answer: Aug 21, 2014 at 21:33 Ismael EL ATIFI
    vector<int> getPermutation(int n, long k) {
        set<int> nums;
        for(int i = 1 ; i <= n ; i++)
            nums.insert(i);
        
        vector<int> result;

        while (nums.size() > 1) {
            long sizeGroup = fact()[nums.size() - 1];
            
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
    set<int> uniqueKs;
    set<int>::iterator kIt;
};

// Heap's algorithm for generating all permutations of 1..n
// adapted from Wikipedia as an iterator
// Allows for infinite iteration
class AllPermsIterator: public PermutationBase {
public:
    AllPermsIterator(int n) : A(n), c(n) {
        setup(n);
    }

    void setup(int n) {
        for (int j = 0; j < n; j++)
            A[j] = j + 1;

        for (int j = 0; j < n; j++)
            c[j] = 0;

        i = 1; 
    }

    vector<int> next() override {
        if (i == A.size())
            setup(A.size());    // allow for infinite iteration

        vector<int> output = A;

        while (c[i] >= i && i < A.size()) {
            c[i] = 0;
            i++;
        }

        if (i < A.size()) {
            if (i % 2 == 0)
                swap(A[0], A[i]);
            else
                swap(A[c[i]], A[i]);
            c[i]++;
            i = 1;
        }

        return output;
    }

private:
    vector<int> A;
    vector<int> c;
    int i;
    vector<int> output;
};

/* 
   Iterator over permutations of 1..n where this sequence is possibly 
   partitioned into m blocks: [1, b1], (b1, b2], .., (bm-1, bm = n].
   Chooses correct iterator for each partition defined by boundaries
   (ie AllPermsIterator or RandomPermsIterator)

*/
class PermutationIterator {
public:
        /*
            If boundaries are specified as {b1, b2, .., bm} then the 
            permutation of 1..bm is partitioned 
            into m blocks: [1, b1], (b1, b2], .., (bm-1, bm]
            
            `count` is then restricted to be no more than max((b[i] - b[i-1])!)

            @param boundaries Right boundaries of NaN blocks (1-based)
        */
    PermutationIterator(long count, vector<int>& boundaries) {
        if (boundaries.size() == 1) {
            setup(count, boundaries[0]);
            return;
        }

        this->boundaries = boundaries ; 
        
        int max_partition_width = 0;
        for (int i = 0; i < boundaries.size(); i++) {
            int partition_width = boundaries[i] - (i == 0 ? 0 : boundaries[i - 1]);
            max_partition_width = max(max_partition_width, partition_width);
        }

        this->count = min(count, PermutationBase::fact()[max_partition_width]); 

        for (int i = 0; i < boundaries.size(); i++) {
            int partition_width = boundaries[i] - (i == 0 ? 0 : boundaries[i - 1]);
            int n = PermutationBase::fact()[partition_width];
            if (count >= n)
                iterators.push_back(new AllPermsIterator(partition_width));
            else {
                iterators.push_back(new RandomPermsIterator(partition_width, count));
            }
        }
    }

        // No boundaries specified, just use 1..n
    PermutationIterator(long count, int n) { setup(count, n); }

        // No count specified, just use n! of 1..n
    PermutationIterator(int n) { setup(PermutationBase::fact()[n], n); }

        // No partitions, just sample 1..n!
    void setup(long count, int n) {
        this->count = count;

        if (count >= PermutationBase::fact()[n])
            iterators.push_back(new AllPermsIterator(n));
        else
            iterators.push_back(new RandomPermsIterator(n, count));
    }

    bool hasNext() const {return count > 0;}

    // Simply concatenate the next permutation from each partition
    // but add boundary[i-1] to each value in partition i
    vector<int> next() {
        vector<int> result;
        for (int i = 0; i < iterators.size(); i++) {
            vector<int> partition_perm = iterators[i]->next();
            for (int v : partition_perm)
                result.push_back((i == 0 ? 0 : boundaries[i-1]) + v);
        }

        count--;

        return result;
    }

    ~PermutationIterator() {
        for (PermutationBase *p : iterators)
            delete p;
    }

private:
    long count;
    vector<int> boundaries; // right boundaries of NaN blocks

    vector<PermutationBase *> iterators;
};

#endif // PERMUTATION_HPP