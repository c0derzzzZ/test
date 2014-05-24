#include <iostream>
#include <string>
#include <exception>
#include <ostream>
#include <iterator>
#include <cctype>
#include <algorithm>
#include <vector>
#include <cctype>
#include <sstream>
#include <deque>
#include <ctime>
#include <thread>
#include <future>
#include <functional>

/*
 * Goal: a minimal big integer class
 *       guidelines from books (exception-safety / allocator / comparator / iterator )
 *       c++11 features
 *       class design
 */


struct Sign final {
  int8_t sign;
};


template<typename T, typename C> class BigNumber;

template<typename T, typename C>
class OperatorIf {
  public:

  // implement assignment version as class member
  // and non-assignment version as non-member using the assignment version
  virtual BigNumber<T,C>& operator += (const BigNumber<T,C>&) = 0;
  virtual BigNumber<T,C>& operator -= (const BigNumber<T,C>&) = 0;
  virtual BigNumber<T,C>& operator *= (const BigNumber<T,C>&) = 0;
//  virtual BigNumber<T,C>& operator /= (const BigNumber<T,C>&) = 0;
//  virtual BigNumber<T,C>& operator %= (const BigNumber<T,C>&) = 0;
//  virtual BigNumber<T,C>& operator ^= (const BigNumber<T,C>&) = 0;
//
//  virtual BigNumber<T,C>& operator ++();    // pre-increment
//
//  // post-increment, return const object to avoid case for e.g. a++++
//  virtual const BigNumber<T,C> operator ++(int);

  // XXX
  // stream operators must be outside the class member
  // std::ostream& operator << (std::ostream& oss, const BigNumber<T,C>& n);
  // std::istream& operator >> (std::istream& oss, BigNumber<T,C>& n);
  //
  // alternative sol:
  virtual std::ostream& print (std::ostream&) const = 0;
  virtual std::istream& read  (std::istream&) = 0;


  virtual bool operator == (const BigNumber<T,C>&) const = 0;
  virtual bool operator != (const BigNumber<T,C>&) const = 0;
  virtual bool operator <  (const BigNumber<T,C>&) const = 0;
//  virtual bool operator >  (const BigNumber<T,C>&) = 0;

  // always make destructor virtual
//  virtual ~OperatorIf();

};

template<typename T>
struct is_bignum {
  static const bool value = false;
};

#if 0
template<typename T, typename C>
struct bignum_traits {
  typedef T digit_type;
  typedef C digit_container;
  typedef
};

template<>
struct is_bignum<T> : public std::true_type {};
#endif

template<typename T=int, typename C=std::vector<T>>
class BigNumber : public OperatorIf<T,C> {
  public:

//  typedef C::iterator       iterator;
//  typedef C::const_iterator const_iterator;

  typedef BigNumber<T,C> self_type;


  // Big 5
  // explicit to avoid implicit conversion and thereby creation of temporaries
  // std::enable_if to allow contructors with only integral types
  // only giving default params doesn't work, give template default as well like int
  template<typename _Tp,
           typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
  /*explicit*/ BigNumber(_Tp value) {
    //std::cout << "is_integral" << std::endl;
    T ten = T(10);
    const _Tp zero = _Tp(0);
    do {
      mDigits.emplace_back(value % ten);
      value /= ten;
    } while (value > zero);

    mSign = 1;
    if (value < zero) mSign = -1;
  }

  BigNumber() : mSign(1), mDigits() {
  //  std::cout << "Default..." << std::endl;
  }

  explicit BigNumber(const std::string& v) { std::cout << "string..." << std::endl;
    auto start = v.rbegin(), end = v.rend();
    if (v.size() > 0 && v[0] == '-') {
      mSign = -1;
      --end;
    } else {
      mSign = 1;
    }

    std::transform(start,
                   end,
                   std::back_inserter(mDigits),
                   [](const char c){ return c - '0'; });
  }
  // copy constructor shouldn't have explicit
  BigNumber (const self_type& rhs) noexcept {
    //std::cout << "copy constructors..." << std::endl;
    std::copy(std::begin(rhs.mDigits), std::end(rhs.mDigits), std::back_inserter(mDigits));
    mSign = rhs.mSign;
  }

  self_type& operator = (const self_type& rhs) noexcept {
    std::cout << "copy assignment.." << std::endl;
    if (this != &rhs) {
      this->mDigits = rhs.mDigits;
      this->mSign = rhs.mSign;
    }
    return *this;
  }
  // TODO: rvalue refrences
  BigNumber (const BigNumber<T,C>&& rhs) noexcept {
    //std::cout << "move construcotr..." << std::endl;
    mSign   = std::move(rhs.mSign);
    mDigits = std::move(rhs.mDigits);
  }

  self_type& operator = (const self_type&& rhs) noexcept {
    //std::cout << "move assignment..." << std::endl;
    if (this != &rhs) {
      this->mSign   = std::move(rhs.mSign);
      this->mDigits = std::move(rhs.mDigits);
    }
    return *this;
  }

  std::string getString () const {
    std::stringstream os;
    //std::cout << "here...." << std::endl;
    if (mSign < 0) os << '-';
    std::copy(mDigits.rbegin(),
              mDigits.rend(),
              std::ostream_iterator<T>(os, ""));
    return os.str();
  }

  std::ostream& print (std::ostream& os) const {
//    std::cout << "here in ostream one.." << std::endl;
    if (mSign < 0) os << '-';
    std::copy(mDigits.rbegin(),
              mDigits.rend(),
              std::ostream_iterator<T>(os, ""));
    return os;
  }

  std::istream& read (std::istream& is) {
    std::string tmp;
    is >> tmp;
    //std::copy(std::istream_iterator<char>(is), std::istream_iterator<char>(), back_inserter(tmp));

    (*this) = self_type(tmp);
    return is;
  }

  size_t size() const {
    return mDigits.size();
  }

  int getInt() {
//    std::stringstream ss;
//    std::cout << "getString=" << getString() << std::endl;
//    ss << getString();
//    int t;
//    ss >> t;
//    if (ss.fail()) std::cerr << "failed" << ss.failbit << std::endl;
    return std::atoi(getString().c_str());
  }

  int getInt() const {
//    std::stringstream ss(getString());
//    std::cout << "getString2 =" << getString() << std::endl;
//    ss << getString();
//    std::cout << "getString22=" << ss.str() << std::endl;
//    std::string t2;
//    ss >> t2;
//    if (ss.fail()) std::cerr << "failed" << ss.failbit << std::endl, ss.clear();
    return std::atoi(getString().c_str());
  }

  bool operator == (const self_type& rhs) const {
    if (rhs.size() != size()) return false;
    if (rhs.mSign != mSign) return false;
    return std::equal(mDigits.begin(), mDigits.end(), rhs.mDigits.begin());
  }

  bool operator != (const self_type& a) const {
    return !operator == (a);
  }

  bool operator <  (const self_type& a) const {
    if (mSign < a.mSign) return true;
    else if (mSign > a.mSign) return false;
    else if (mSign > 0) {
      if (size() < a.size()) return true;
      if (size() > a.size()) return false;
      //std::cout << "size..." << size() << ", " << a.size() << std::endl;
      return std::lexicographical_compare(mDigits.rbegin(), mDigits.rend(), a.mDigits.rbegin(), a.mDigits.rend());
    } else {
      if (size() > a.size()) return true;
      if (size() < a.size()) return false;
      return !std::lexicographical_compare(mDigits.rbegin(), mDigits.rend(), a.mDigits.rbegin(), a.mDigits.rend());
    }
  }

  self_type& operator -() {
    mSign = -mSign;
    return *this;
  }

  self_type& operator += (const self_type& rhs) {
    if (mSign == rhs.mSign) return sumSameType(rhs);
//    else if ((*this) < rhs) return rhs - (*this);
//    else return (*this) - rhs;
    return *this;
  }

  self_type& operator -= (const BigNumber<T,C>& rhs) {
//    std::cout << "operator -" << std::endl;
    if (mSign != rhs.mSign) return sumSameType(rhs);
    return subtract(rhs);
  }

  self_type& operator *= (const self_type& rhs) {
    self_type tmp(*this);
    return (*this) = karatsuba(tmp, rhs);
  }

  virtual ~BigNumber() noexcept {}

  // friendship doesn't inherit
  template<typename TT, typename CC>
  friend std::ostream& operator << (std::ostream& oss, const BigNumber<TT,CC>& n);

  template<typename TT, typename CC>
  friend std::istream& operator >> (std::istream& is, BigNumber<TT,CC>& n);

//  template<typename TT, typename CC>
//  friend const BigNumber<T,C> operator + (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);
//
//  template<typename TT, typename CC>
//  friend const BigNumber<TT,CC> operator - (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);
//
//  template<typename TT, typename CC>
//  friend const BigNumber<TT,CC> operator * (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);
//
//  template<typename TT, typename CC>
//  friend const BigNumber<TT,CC> operator / (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);
//
//  template<typename TT, typename CC>
//  friend const BigNumber<TT,CC> operator % (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);
//
//  template<typename TT, typename CC>
//  friend const BigNumber<TT,CC> operator ^ (const BigNumber<TT,CC>&, const BigNumber<TT,CC>&);

  protected:

  self_type& subtract(const self_type& rhs) {
    T carry = T(0);
    C tmp;
    auto myNumDigits = mDigits.size(), rhsNumDigits = rhs.mDigits.size();
    auto n = std::max(myNumDigits, rhsNumDigits), i = 0u;

    auto digitSubtractor = [&](T a, T b) {
      auto tmpA = a+carry;
      if (b <= tmpA) {
        carry = T(0);
        return tmpA-b;
      } else {
        carry = T(-1);
        return T(10)+tmpA-b;
      }
    };

    if ((*this) < rhs) {
      for (auto it = mDigits.begin(), jt = rhs.mDigits.begin(); i < n; ++i) {
        T a = T(0), b = T(0);
        if (it != mDigits.end()) a = *it, ++it;
        if (jt != rhs.mDigits.end()) b = *jt, ++jt;
        tmp.emplace_back(digitSubtractor(b, a));
      }
      mSign = -1;
      //if (carry < T(0)) tmp.emplace_back(carry);
    } else {
      for (auto it = mDigits.begin(), jt = rhs.mDigits.begin(); i < n; ++i) {
        T a = T(0), b = T(0);
        if (it != mDigits.end()) a = *it, ++it;
        if (jt != rhs.mDigits.end()) b = *jt, ++jt;
        tmp.emplace_back(digitSubtractor(a, b));
      }
      mSign = +1;
      //if (carry < T(0)) tmp.emplace_back(carry);
    }

    mDigits.swap(tmp);
    trim();
    return (*this);
  }

  self_type& sumSameType(const self_type& rhs) {
    T carry = T(0);
    C tmp;
    auto myNumDigits = mDigits.size(), rhsNumDigits = rhs.mDigits.size();
    auto n = std::max(myNumDigits, rhsNumDigits), i = 0u;

    auto digitAdder = [&](T a, T b) {
      T sum = a+b+carry;

      //std::cout << a << "+" << b << "+" << carry << "=" << sum << std::endl;
      if (sum >= T(10)) carry = T(1);
      else carry = T(0);
      return sum % T(10);
    };


    for (auto it = mDigits.begin(), jt = rhs.mDigits.begin(); i < n; ++i) {
      T a = T(0), b = T(0);
      if (it != mDigits.end()) a = *it, ++it;
      if (jt != rhs.mDigits.end()) b = *jt, ++jt;
      tmp.emplace_back(digitAdder(a, b));
    }
    if (carry > T(0)) tmp.emplace_back(carry);

    //std::cout << "tmp = "; std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<T>(std::cout, "")); std::cout << std::endl;
    mDigits.swap(tmp);
    return (*this);
  }



  self_type karatsuba(const self_type& a, const self_type& b) {
//    std::cout << "karatsuba: " << a << "*" << b << std::endl;
    size_t alen = a.size(), blen = b.size();
    size_t m = std::max(alen, blen);

    if (m < 8 && ((a < self_type::mTen) || (b < self_type::mTen))) {
      //std::cout << "ten = " << self_type::mTen << std::endl;
      int v1 = a.getInt();
      int v2 = b.getInt();
//      std::cout << "base case: " << a << "*" << b << "=" << (v1*v2) << std::endl << std::endl;
      return (v1*v2);
    }

    size_t m2 = m/2;
    C cl1(a.mDigits.begin(),    a.mDigits.begin()+m2),
              ch1(a.mDigits.begin()+m2, a.mDigits.end()),
              cl2(b.mDigits.begin(),    b.mDigits.begin()+m2),
              ch2(b.mDigits.begin()+m2, b.mDigits.end());

    self_type l1,h1,l2,h2;
    l1.mDigits.swap(cl1);
    h1.mDigits.swap(ch1);
    l2.mDigits.swap(cl2);
    h2.mDigits.swap(ch2);

    const typename C::value_type zero = typename C::value_type(0);
//    h1.mDigits.resize(m, zero);
//    l1.mDigits.resize(m, zero);
//    h2.mDigits.resize(m, zero);
//    l2.mDigits.resize(m, zero);


//    for (size_t i = 0; i < m2; ++i) {
//      if (i < alen) l1.mDigits.push_back(a.mDigits[i]);
//      if (i < blen) l2.mDigits.push_back(b.mDigits[i]);
//    }
//
//    for (size_t i = m2; i < m; ++i) {
//      if (i < alen) h1.mDigits.push_back(a.mDigits[i]);
//      if (i < blen) h2.mDigits.push_back(b.mDigits[i]);
//    }

//    std::cout << "l1 = " << l1 << std::endl;
//    std::cout << "h1 = " << h1 << std::endl;
//    std::cout << "l2 = " << l2 << std::endl;
//    std::cout << "h2 = " << h2 << std::endl;

//    self_type ta = l1+h1, tb = l2+h2;
//    std::packaged_task<self_type()> p0(std::bind(&self_type::karatsuba, this, std::cref(l1), std::cref(l2)));
//    std::packaged_task<self_type()> p2(std::bind(&self_type::karatsuba, this, std::cref(h1), std::cref(h2)));
//    std::packaged_task<self_type()> p1(std::bind(&self_type::karatsuba, this, std::cref(ta), std::cref(tb)));
//
//    std::future<self_type> r0 = p0.get_future();
//    std::future<self_type> r2 = p2.get_future();
//    std::future<self_type> r1 = p1.get_future();
//
//    p0();
//    p1();
//    p2();
//
//    self_type z0 = r0.get();
//    self_type z2 = r2.get();
//    self_type z1 = r1.get();

        self_type z0 = karatsuba(l1, l2);
        self_type z2 = karatsuba(h1, h2);
        self_type z1 = karatsuba(l1+h1, l2+h2);

    // z2*10^(2*m2))+((z1-z2-z0)*10^(m2))+(z0)
    self_type y = z1 - (z2 + z0);

//    std::cout << "z0 = " << z0 << std::endl;
//    std::cout << "z1 = " << z1 << std::endl;
//    std::cout << "z2 = " << z2 << std::endl;
//    std::cout << "y  = " << y << std::endl;
//    std::cout << "m = " << m << ", m2 = " << m2 << std::endl;
    C z2Tmp(2*m2, zero);
    //std::fill_n(begin(z2Tmp), 2*m2, zero);
    std::copy(begin(z2.mDigits), end(z2.mDigits), std::back_inserter(z2Tmp));
    z2.mDigits.swap(z2Tmp);
//    std::cout << "new z2 = " << z2 << std::endl;

    C yTmp(m2, zero);
    //std::fill_n(begin(yTmp), m2, zero);    // xxx: fill_n need space alreay available
    std::copy(begin(y.mDigits), end(y.mDigits), std::back_inserter(yTmp));
    y.mDigits.swap(yTmp);

    //y.trim();
    //z2.trim();
    self_type ret = z2 + y + z0;
    z2.trim();
//    std::cout << "karatsuba: " << a << "*" << b << " = " << ret << std::endl << std::endl;
    return ret;
  }

  void trim() {
    auto it = std::find_if(mDigits.rbegin(), mDigits.rend(), [](T& v) { return (v != T(0));});
//    std::cout << "digit = " << *it << std::endl;
    auto dis = std::distance(it, mDigits.rend());
    auto s = mDigits.begin();
    std::advance(s, dis);
//    std::cout << "tmp = "; std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<T>(std::cout, "")); std::cout << std::endl;
//    std::cout << "dis = " << dis << std::endl;
    mDigits.erase(s, mDigits.end());
    if (mDigits.empty()) { mDigits.push_back(T(0)); }
//    std::cout << "tmp = "; std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<T>(std::cout, "")); std::cout << std::endl;
  }

  private:

  C   mDigits;   // least to most significant
  int mSign;
  static const self_type mZero, mTen;
};

template<typename TT, typename CC>
std::ostream& operator << (std::ostream& os, const BigNumber<TT,CC>& n) {
  return n.print(os);
}

template<typename T, typename C>
std::istream& operator >> (std::istream& is, BigNumber<T,C>& n) {
  return n.read(is);
}

template<typename _Tp,
  typename = typename std::enable_if<std::is_integral<_Tp>::value || is_bignum<_Tp>::value>>
const _Tp operator + (const _Tp& a, const _Tp& b) {
  _Tp tmp(a);
  tmp += b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator + (const _Tp& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp += b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator + (const BigNumber<T,C>& b, const _Tp& a) {
  BigNumber<T,C> tmp(a);
  tmp += b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator + (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp += b;
  return tmp;
}

template<typename _Tp,
  typename = typename std::enable_if<std::is_integral<_Tp>::value || is_bignum<_Tp>::value>>
const _Tp operator - (const _Tp& a, const _Tp& b) {
  _Tp tmp(a);
  tmp -= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator - (const _Tp& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp -= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator - (const BigNumber<T,C>& b, const _Tp& a) {
  BigNumber<T,C> tmp(a);
  tmp -= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator - (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp -= b;
  return tmp;
}

// multiplicaiton
template<typename _Tp,
  typename = typename std::enable_if<std::is_integral<_Tp>::value || is_bignum<_Tp>::value>>
const _Tp operator * (const _Tp& a, const _Tp& b) {
  _Tp tmp(a);
  tmp *= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator * (const _Tp& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp *= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator * (const BigNumber<T,C>& b, const _Tp& a) {
  BigNumber<T,C> tmp(a);
  tmp *= b;
  return tmp;
}

template<typename T, typename C, typename _Tp,
typename = typename std::enable_if<std::is_integral<_Tp>::value, _Tp>::type>
const BigNumber<T,C> operator * (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
  BigNumber<T,C> tmp(a);
  tmp *= b;
  return tmp;
}


//template<typename T, typename C>
//friend const BigNumber<T,C> operator - (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
//  BigNumber<T,C> tmp(a);
//  tmp -= b;
//  return tmp;
//}
//
//template<typename T, typename C>
//friend const BigNumber<T,C> operator * (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
//  BigNumber<T,C> tmp(a);
//  tmp *= b;
//  return tmp;
//}
//
//template<typename T, typename C>
//friend const BigNumber<T,C> operator / (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
//  BigNumber<T,C> tmp(a);
//  tmp /= b;
//  return tmp;
//}
//
//template<typename T, typename C>
//friend const BigNumber<T,C> operator % (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
//  BigNumber<T,C> tmp(a);
//  tmp %= b;
//  return tmp;
//}
//
//template<typename T, typename C>
//friend const BigNumber<T,C> operator ^ (const BigNumber<T,C>& a, const BigNumber<T,C>& b) {
//  BigNumber<T,C> tmp(a);
//  tmp ^= b;
//  return tmp;
//}

template<typename T, typename C>
const BigNumber<T,C> BigNumber<T,C>::mZero(0);

template<typename T, typename C>
const BigNumber<T,C> BigNumber<T,C>::mTen(10);

typedef BigNumber<int, std::vector<int>> BigInt;

int main() {
//  BigInt a(100u); // is_integral contors
//  std::cout << std::endl;
//  BigInt b = a;   // copy contors
//  std::vector<BigInt> vec;
//  std::cout << std::endl;
//  vec.emplace_back(10000);  // in-place construction (is_integral)
//  std::cout << std::endl;
//  vec.emplace_back(BigInt("1"));
//  std::cout << std::endl;
//  vec.push_back(a);
//  std::cout << std::endl;
//  vec.push_back(b);
//  std::cout << std::endl;
//  std::copy(begin(vec), end(vec), std::ostream_iterator<BigInt>(std::cout, "\n"));
  std::clock_t start = clock();
  BigInt x("999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555");
  //BigInt y("10000000000000000000000000000000000000000000000000000000000000000000000000000");
  BigInt y("9999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555999999999999999999999999999999999999999999999999999999999999999999999999999991222222222222222222223333333333333333333333399999999999999999999999999999999999999999999999888888888888888888888888888883333333333333333333333333333333333333333333333333333333333333339999999999999991111111111133333333333335555555555555555555555555555555555555555555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999999122222222222222222222333333333333333333333339999999999999999999999999999999999999999999999988888888888888888888888888888333333333333333333333333333333333333333333333333333333333333333999999999999999111111111113333333333333555555555555555555555555555555555555555555555555555555555555555559999999999999999999999999999999999999999999999999999999999999999999999999999912222222222222222222233333333333333333333333999999999999999999999999999999999999999999999998888888888888888888888888888833333333333333333333333333333333333333333333333333333333333333399999999999999911111111111333333333333355555555555555555555555555555555555555555555555555555555555555555");
//  std::cout << x << std::endl
//            << y << std::endl;
//  std::cout << x-y << std::endl;
  std::cout << "ans = " << x*y << std::endl;
  //std::cout << x << "+" << y << "=" << -(x+y+1) << std::endl;
  std::clock_t end = clock();
  fprintf(stderr, "time=%.3lfsec\n", 0.001 * (clock() - start));
#if 0
  while (true) {
    std::cin >> a >> b;
    std::cout << a << ", " << b << std::endl;
    if (b < a) std::cout << "b < a" << std::endl;
  }
#endif
}
