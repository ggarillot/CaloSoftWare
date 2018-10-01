#include <iostream>
#include <algorithm>

namespace algorithm{
  
  template<class T>
    bool comparePairValues(const T & x, const T & y)
    {
      return x.second < y.second;
    }

  template<class T>
    typename T::const_iterator map_max_element(const T &A)
    {
      return std::max_element(A.begin(), A.end(), comparePairValues<typename T::value_type>);
    }

  template<class T>
    typename T::const_iterator map_min_element(const T &A)
    {
      return std::min_element(A.begin(), A.end(), comparePairValues<typename T::value_type>);
    }

}
