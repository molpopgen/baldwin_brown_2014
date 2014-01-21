#ifndef __INITMS_HPP__
#define __INITMS_HPP__
#include <isbinary.hpp>
#include <forward_types.hpp>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include <limits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <iostream>
#include <Sequence/SimData.hpp>

template< typename gamete_type,
	  typename vector_allocator_type,
	  typename mutation_type,
	  typename list_allocator_type,
	  template<typename,typename> class vector_type,
	  template<typename,typename> class list_type >
void init_with_ms( gsl_rng * r,
		   Sequence::SimData & d,
		   vector_type<gamete_type,vector_allocator_type> * gametes,
		   list_type<mutation_type,list_allocator_type> * mutations,
		   const unsigned & max_chroms = std::numeric_limits<unsigned>::max())
{
  typedef list_type<mutation_type,list_allocator_type> MLIST_TYPE;
  /*
    add the mutations to the list.
    We create a lookup table of position -> pointer to mutations,
    to keep the search time logarithmic rather than linear 
    later on.   This makes a HUGE speed difference.
  */
  std::map<double,typename  MLIST_TYPE::iterator> mutlookup;

  unsigned site=0;
  for( Sequence::SimData::pos_iterator i = d.pbegin() ; 
       i != d.pend() ; ++i,++site )
    {
      //is it a seg site in the first max_chroms?
      unsigned c=0;
      for(unsigned ind=0;ind<std::min(max_chroms,unsigned(d.size()));++ind)
	{
	  if((d)[ind][site]=='1')++c;
	}
      if(i != d.pend()-2)
	{
	  if(*(i+1) == *(i+2))
	    {
	      //new position value is 1/2 way between i and i+1.
	      //std::cerr << *i << ' ' << *(i+1) << ' ' <<*(i+2) << "->";
	      *(i+1) = ( (*i) + ( (*(i+1)-*i)/2. ) );
	      //std::cerr << *i << ' ' << *(i+1) << ' ' << *(i+2) << "\n";
	    }
	}
      //mutations->push_back( mutation_type(*i,0,0) );
      if(c)
	{
	  typename MLIST_TYPE::iterator itr = mutations->insert(mutations->end(),
								mutation_type(*i,0,0));
	  assert( mutlookup.find(*i) == mutlookup.end() );
	  mutlookup.insert(std::make_pair(*i,itr));
	}
    }
  //Figure out the unique set of gametes, and how frequent they are
  std::set<std::string> ugams(d.begin(),std::min(d.begin()+max_chroms+1,d.end()));
  for(std::set<std::string>::const_iterator i = ugams.begin(); i != ugams.end() ; ++i )
    {
      unsigned c = std::count(d.begin(),std::min(d.begin()+max_chroms,d.end()),*i);
      //now, create this gamete and fill its neutral mutations list
      gamete_type g(c);
      for( std::string::const_iterator j = i->begin() ; j != i->end() ; ++j )
	{
	  if( *j == '1' )//is a derived mutation
	    {
	      double pos = *(d.pbegin() + (j-i->begin()));
	      //find the mutation in the mutations list
	      if( mutlookup.find(pos) != mutlookup.end() )
		{
		  mutlookup[pos]->n += c;
		  g.mutations.push_back(mutlookup[pos]);
		}
	    }
	}
      gametes->push_back(g);
    }
  ugams.erase(ugams.begin(),ugams.end());
}

template< typename streamtype,
	  typename gamete_type,
	  typename vector_allocator_type,
	  typename mutation_type,
	  typename list_allocator_type,
	  template<typename,typename> class vector_type,
	  template<typename,typename> class list_type >
void init_with_ms_from_stream( gsl_rng * r,
			       streamtype & in,
			       vector_type<gamete_type,vector_allocator_type> * gametes,
			       list_type<mutation_type,list_allocator_type> * mutations,
			       const unsigned & max_chroms = std::numeric_limits<unsigned>::max())
{
  typedef list_type<mutation_type,list_allocator_type> MLIST_TYPE;
  Sequence::SimData d;
  in >> d >> std::ws;
  init_with_ms( r,d,gametes,mutations,max_chroms );
  /*
    add the mutations to the list.
    We create a lookup table of position -> pointer to mutations,
    to keep the search time logarithmic rather than linear 
    later on.   This makes a HUGE speed difference.
  */
  /*
  std::map<double,typename  MLIST_TYPE::iterator> mutlookup;

  unsigned site=0;
  for( Sequence::SimData::pos_iterator i = d->pbegin() ; 
       i != d->pend() ; ++i,++site )
    {
      //is it a seg site in the first max_chroms?
      unsigned c=0;
      for(unsigned ind=0;ind<std::min(max_chroms,unsigned(d->size()));++ind)
	{
	  if((*d)[ind][site]=='1')++c;
	}
      if(i != d->pend()-2)
	{
	  if(*(i+1) == *(i+2))
	    {
	      //new position value is 1/2 way between i and i+1.
	      //std::cerr << *i << ' ' << *(i+1) << ' ' <<*(i+2) << "->";
	      *(i+1) = ( (*i) + ( (*(i+1)-*i)/2. ) );
	      //std::cerr << *i << ' ' << *(i+1) << ' ' << *(i+2) << "\n";
	    }
	}
      //mutations->push_back( mutation_type(*i,0,0) );
      if(c)
	{
	  typename MLIST_TYPE::iterator itr = mutations->insert(mutations->end(),
								mutation_type(*i,0,0));
	  assert( mutlookup.find(*i) == mutlookup.end() );
	  mutlookup.insert(std::make_pair(*i,itr));
	}
    }
  //Figure out the unique set of gametes, and how frequent they are
  std::set<std::string> ugams(d->begin(),std::min(d->begin()+max_chroms+1,d->end()));
  for(std::set<std::string>::const_iterator i = ugams.begin(); i != ugams.end() ; ++i )
    {
      unsigned c = std::count(d->begin(),std::min(d->begin()+max_chroms,d->end()),*i);
      //now, create this gamete and fill its neutral mutations list
      gamete_type g(c);
      for( std::string::const_iterator j = i->begin() ; j != i->end() ; ++j )
	{
	  if( *j == '1' )//is a derived mutation
	    {
	      double pos = *(d->pbegin() + (j-i->begin()));
	      //find the mutation in the mutations list
	      if( mutlookup.find(pos) != mutlookup.end() )
		{
		  mutlookup[pos]->n += c;
		  g.mutations.push_back(mutlookup[pos]);
		}
	    }
	}
      gametes->push_back(g);
    }
  ugams.erase(ugams.begin(),ugams.end());
  delete d;
  */
}

#endif
