#ifndef dout_tools_h
#define dout_tools_h

#include <iostream>
using namespace std;

bool debug=false;

void dout() { std::cout << std::endl; }
template <typename Head, typename... Tail>
void dout(Head H, Tail... T) {
    std::cout << H << ' ';
    dout(T...);
}

void doutnonl() {}
template <typename Head, typename... Tail>
void doutnonl(Head H, Tail... T) {
    std::cout << H << ' ';
    doutnonl(T...);
}


void ddoutnonl() {}
template <typename Head, typename... Tail>
void ddoutnonl(Head H, Tail... T) {
    if (debug) std::cout << H << ' ';
    ddoutnonl(T...);
}

void ddout() { if (debug) std::cout << std::endl; }
template <typename Head, typename... Tail>
void ddout(Head H, Tail... T) {
    if (debug) dout(H, T...);
}

template<typename T>
void printVector(const T& t) {
    if (t.size() > 0)
    {
        std::copy(t.cbegin(), t.cend() - 1, std::ostream_iterator<typename T::value_type>(std::cout, ", "));
        dout(t.back());
    }
    else dout("{}");
}
template<typename T>
void printVectorInVector(const T& t) {
    std::for_each(t.cbegin(), t.cend(), printVector<typename T::value_type>);
}

template<typename T>
void dprintVector(const T& t) {
    if (debug)
    {
        if (t.size() > 0)
        {
          std::copy(t.cbegin(), t.cend() -1, std::ostream_iterator<typename T::value_type>(std::cout, ", "));
          dout(t.back());
        }
        else dout("{}");
    }
}
template<typename T>
void dprintVectorInVector(const T& t) {
    if (debug) std::for_each(t.cbegin(), t.cend(), dprintVector<typename T::value_type>);
}

template <class T>
bool contains(std::vector<T> const &v, T const &x) {
    return ! (v.empty() && std::find(v.begin(), v.end(), x) == v.end());
}

#endif
