//  -*- C++ -*- 
#ifndef _HAPLOID_FUNCTIONS_TCC_
#define _HAPLOID_FUNCTIONS_TCC_
#include <utility>
#include <functional>
#include <vector>
#include <util.hpp>
#ifndef NO_KAHAN_SUMS
#include <kahan.hpp>
#endif

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename haploid_fitness_function>
double sample_haploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & N, const haploid_fitness_function & ff)
{
#ifdef NO_KAHAN_SUMS
  std::vector<double> expected_gamete_freqs(gametes->size(),0.);
  double wbar = 0.;
#else
  std::vector<kahan_sum<double> > expected_gamete_freqs(gametes->size(),kahan_sum<double>());
  kahan_sum<double> wbar;
#endif
  for(unsigned i=0;i<gametes->size();++i)
    {
      double w = ff( *(gametes->begin()+i) );
      expected_gamete_freqs[i] += (w*double((gametes->begin()+i)->n)/double(N));
      wbar += (w*double((gametes->begin()+i)->n)/double(N));
    }
	
#ifdef NO_KAHAN_SUMS
  std::transform(expected_gamete_freqs.begin(),expected_gamete_freqs.end(),expected_gamete_freqs.begin(),std::bind2nd(std::divides<double>(),wbar));
#else
  std::transform(expected_gamete_freqs.begin(),expected_gamete_freqs.end(),expected_gamete_freqs.begin(),std::bind2nd(std::divides<kahan_sum<double> >(),wbar));
#endif
  multinomial_sample_gametes(r, gametes->begin(), gametes->size(), &expected_gamete_freqs[0],N);

  update_gamete_list(gametes,N);
  return wbar;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename haploid_fitness_function>
double sample_haploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, const unsigned & Ncurr, const unsigned & Nnext, const haploid_fitness_function & ff)
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
#ifdef NO_KAHAN_SUMS
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
#else
  kahan_sum<double> wbar;
  std::vector<kahan_sum<double> > expected_gamete_frequencies(gametes->size(),kahan_sum<double>());
#endif
    for(unsigned i=0;i<gametes->size();++i)
    {
      double w = ff( *(gametes->begin()+i) );
      expected_gamete_frequencies[i] += (w*double((gametes->begin()+i)->n)/double(Ncurr));
      wbar += (w*double((gametes->begin()+i)->n)/double(Ncurr));
    }
  //normalize expected frequqncies
#ifdef NO_KAHAN_SUMS
  transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
#else
  transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<kahan_sum<double> >(),wbar));						
#endif
  multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],Nnext);
  update_gamete_list(gametes,Nnext);
  return wbar;
}
#endif /* _HAPLOID_FUNCTIONS_TCC_ */
