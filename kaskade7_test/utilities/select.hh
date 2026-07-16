#ifndef SELECT_HH
#define SELECT_HH

/// Deprecated (use std::conditional_t instead). Will be removed after 2019-12-31
template <bool decide, class T1, class T2>
struct Select{
  typedef T1 type;
};

template <class T1, class T2>
struct Select<false,T1,T2>{
  typedef T2 type;
};

#endif
