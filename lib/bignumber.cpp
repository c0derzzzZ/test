
// Design patterns

template<typename _Traits>
class Number {
  public:

  typedef _Traits                           traits_type;
  typedef typename _Traits::word_type       word_type;
  typedef typename _Traits::container_type  container_type;

  typedef Number<traits_type>   self_type;

  virtual self_type& operator += (const self_type& rhs);

  private:


};

template<typename C>
struct Number_traits {};

template<>
struct Number_traits<Number> {

};


// strategy pattern
class Operation {
  virtual const Number execute() const;
  virtual ~Operation() noexcept;
};

class UnaryOperation : public Operation {
  public:

};

class BinaryOperation : public Operation {

};

class Multiplication : public BinaryOperation {
  public:

  virtual ~Multiplication() noexcept;
};

