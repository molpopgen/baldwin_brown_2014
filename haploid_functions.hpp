#ifndef _HAPLOID_FUNCTIONS_HPP_
#define _HAPLOID_FUNCTIONS_HPP_

#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <forward_types.hpp>
#include <util.hpp>

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename haploid_fitness_function>
double sample_haploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, const unsigned & N, const haploid_fitness_function & ff);

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename haploid_fitness_function>
double sample_haploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, const unsigned & Ncurr, const unsigned & Nnext, const haploid_fitness_function & ff);

#endif /* _HAPLOID_FUNCTIONS_HPP_ */
#include <haploid_functions.tcc>
