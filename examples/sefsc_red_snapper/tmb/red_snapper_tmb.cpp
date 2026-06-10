// Placeholder TMB reference implementation for the SEFSC red-snapper-style example.
//
// The next milestone should implement the same likelihood and derived quantities
// as the Quadra model so objective values, estimates, random effects, and
// uncertainty outputs can be compared side by side.

template<class Type>
Type objective_function<Type>::operator()() {
  return Type(0.0);
}
