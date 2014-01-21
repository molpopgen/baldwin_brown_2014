#include <diploid.hh>
#include <initms.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <isbinary.hpp>
#include <sstream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <Sequence/PolyTableFunctions.hpp>

using namespace std;
using namespace Sequence;
using namespace boost::iostreams;

struct breakpoint : public mutation_base
{
  /*
    A gamete will be a list of these breakpoints stored in its neutral mutation vector.
    When sorted by breakpoint position, a gamete of of founder type hap from the position of 
    the previous breakpoint until the position of the current breakpoint.
   */
  mutable unsigned hap; 
  double s,h;
  breakpoint( const double & position, 
	      const unsigned & count,const unsigned & original_type,
	      const double & sel,const double & dom,
	      const bool & n=true) 
    : mutation_base(position,count,n),hap(original_type),s(sel),h(dom)
  {
  }
  bool operator==(const breakpoint & rhs) const
  {
    return( fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	    this->hap == rhs.hap && this->s == rhs.s );
  }
};

typedef breakpoint mtype;
typedef gamete_base<breakpoint> gtype;

#ifndef NDEBUG
bool check_sum(const std::vector<gtype> & gametes, const unsigned & twoN)
{
  unsigned check=0;
  for(std::vector<gtype>::const_iterator i=gametes.begin();i<gametes.end();++i)
    {
      check+=i->n;
    }
  return (check == twoN);
}
#endif


struct find_first : public std::binary_function<std::pair<double,std::string>,double,bool >
{
  inline bool operator()( const std::pair<double,std::string> & p,const double & d)const
  {
    return d==p.first;
  }
};

struct find_first_g : public std::binary_function<std::pair<double,std::string>,double,bool >
{
  inline bool operator()( const std::pair<double,std::string> & p,const double & d)const
  {
    return p.first > d;
  }
};

struct find_first_le : public std::binary_function<std::pair<double,std::string>,double,bool >
{
  inline bool operator()( const std::pair<double,std::string> & p,const double & d)const
  {
    return p.first <= d;
  }
};

struct find_first2 : public std::binary_function<std::pair<double, vector<double> >,double,bool >
{
  inline bool operator()( const pair<double,vector<double> > & lhs, const double & d) const
  {
    return d==lhs.first;
  }
};

struct sort_on_first : public std::binary_function< std::pair<double,std::string>,
						    std::pair<double,std::string>,
						    bool >
{
  inline bool operator()( const std::pair<double,std::string> & lhs,
			  const std::pair<double,std::string> & rhs) const
  {
    return lhs.first < rhs.first;
  }
};

#ifndef NDEBUG
bool mutations_are_sorted( const vector<gtype> & gametes )
{
  for(unsigned i=0;i<gametes.size();++i)
    {
      unsigned j=0,k=1;
      for( ; k < gametes[i].mutations.size() ; ++j,++k )
	{
	  if( gametes[i].mutations[k]->pos < gametes[i].mutations[j]->pos )
	    {
	      return false;
	    }
	}
    }
  return true;
}

bool nmuts_not_empty( const vector<gtype> & gametes )
{
  //this is an incorrect check--breakpoints can fix!!!
  return true;
  for(unsigned i=0;i<gametes.size();++i)
    {
      if( gametes[i].mutations.empty() )
	{
	  cerr << ' ' << gametes[i].n << ' ' << gametes[i].smutations.size() << '\n';
	  return false;
	}
    }
  return true;
}
#endif

struct mutation_le_pos
//finds first breakpoint >= position
{
  typedef bool result_type;
  template<typename mutation_type>
  inline bool operator()(const mutation_type & m, const double & d) const
  {
    return( d >= m.pos );
  }
};

struct mutless : public binary_function<mtype,mtype,bool>
{
  inline bool operator()(const mtype &l,const mtype &r)const
  {
    return l.pos < r.pos;
  }
};


// template< typename gamete_type,
// 	  typename vector_type_allocator,
// 	  template <typename,typename> class vector_type >
// void update_gamete_list2( vector_type<gamete_type,vector_type_allocator > * gametes, 
// 			  const unsigned & N)
struct update_gamete_list2
{
typedef void result_type;
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template <typename,typename> class vector_type >
  void operator()( vector_type<gamete_type,vector_type_allocator > * gametes, 
			    const unsigned & N) const
  {
  int g = gametes->size()-1;
  typename vector_type<gamete_type,vector_type_allocator >::iterator gbeg;
  while(g>=0)
    {
      gbeg=gametes->begin()+g;
      assert(gbeg->n <= N);
      if( gbeg->n == 0 )
	{
	  gametes->erase(gbeg);
	}
      else
	{
	  int j = gbeg->mutations.size()-1;
	  while(j >= 0)
	    {
	      assert( gbeg->mutations[j]->n <= N );
	      //if(gbeg->mutations[j]->n == N|| gbeg->mutations[j]->n == 0)
	      if(gbeg->mutations[j]->n == 0)
		{
		  gbeg->mutations.erase(gbeg->mutations.begin()+j);
		}
	      --j;
	    }
	  j = gbeg->smutations.size()-1;
	  while(j >= 0)
	    {
	      assert( gbeg->smutations[j]->n <= N );
	      //if(gbeg->smutations[j]->n == N|| gbeg->smutations[j]->n == 0)
	      if(gbeg->smutations[j]->n == 0)
		{
		  gbeg->smutations.erase(gbeg->smutations.begin()+j);
		}
	      --j;
	    }
	}
      --g;
    }
}
};

template<typename mutationtype,
	 //typename vector_type_allocator1,
	 //typename vector_type_allocator2,
	 typename list_type_allocator,
	 //template <typename,typename> class vector_type,
	 template <typename,typename> class list_type >
void remove_lost( list_type<mutationtype,list_type_allocator> * mutations, 
		  const unsigned & generation,const unsigned & N)
{
  BOOST_STATIC_ASSERT( ( boost::is_base_and_derived<mutation_base,mutationtype>::value) );
  typename list_type<mutationtype,list_type_allocator>::iterator i = mutations->begin();
  
  while(i != mutations->end())
    {
      assert(i->n <= N);			
      i->checked=false;
      // if(i->n==N )
      // 	{
      // 	  fixations->push_back(*i);
      // 	  fixation_times->push_back(generation);
      // 	}
      if( i->n == 0 ) //|| i->n == N )
	{
	  mutations->erase(i);
	  i=mutations->begin();
	}
      else
	{
	  ++i;
	}
    }
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function,
	  typename updater>
double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		       const unsigned & twoN, const diploid_fitness_function & ff,const updater & u,const double & f = 0.)
  //hermaphroditic diploids
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
  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],gametes,twoN,ff,f);
  //normalize expected frequqncies
#ifdef NO_KAHAN_SUMS
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
#else
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<kahan_sum<double> >(),wbar));						
#endif
  multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],twoN);
    
  //update_gamete_list2(gametes,twoN);
  u(gametes,twoN);
  return wbar;
}

template< typename gamete_type,
	  typename mutation_type,
	  typename vector_type_allocator,
 	  template<typename,typename> class vector_type,
	  typename map_function>
unsigned recombine_breakpoints(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes,
			       list<mutation_type> * mutations,
			       const unsigned & twoN, const double & littler,const map_function & mf,
			       const double & f = 0)
/*
  Note: assuming only 1 selected site @ the moment.  Otherwise, we don't make an guarantees
*/
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

  typedef typename gamete_type::mutation_container mcont;
  typedef typename mcont::iterator mcont_iterator;
  typedef typename vector_type<gamete_type, vector_type_allocator>::iterator vtype_iterator;
  //"individual-based recombination"
  //1. determine # recombinants in whole pop
  //2. assign them randomly to pairs of gametes chosen from the gametes list
  //3. if they are the same gamete, ignore, else do recombination thing.

  //calc freq of homozygotes
  
  double fAA=0.;
  if(gametes->size() < twoN)
    {
      const double twoNsq = double(twoN)*double(twoN);
      for(vtype_iterator i=gametes->begin() ;i<gametes->end();++i)
	{
	  fAA += ( (double(i->n)*double(i->n)*(1.-f))/twoNsq + (double(i->n)/twoN)*f );
	}
    }
  unsigned NRECS = unsigned(gsl_ran_poisson(r,(1.-fAA)*double(twoN)*littler));
  unsigned doublecheck=NRECS;
#ifndef NDEBUG
  unsigned NRECS_DONE = 0;
#endif
  vtype_iterator ibeg,jbeg;
  const unsigned ncurrent_classes = gametes->size();
  while(NRECS > 0)
    {
      pick_gametes(r,gametes,&ibeg,&jbeg, ncurrent_classes);
      unsigned nm1=ibeg->mutations.size()+ibeg->smutations.size(),
	nm2=jbeg->mutations.size()+jbeg->smutations.size();
      //if one gametes carries 0 mutations, and the other carries 1,
      //it is impossible for recombination to generate any new types,
      //so we skip those cases...
      if(!(std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1))
	{
	  ibeg->n--;
	  jbeg->n--;
	      
	  gamete_type new_gamete1(1, mcont(),mcont()),new_gamete2(new_gamete1);
	  double pos = mf();
	  //make sure this position does not exist anywhere ("infinitely-many recombination breakpoints")
	  while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) 
		 != mutations->end() )
	    {
	      pos = mf();
	    }
	  //pointer arithmetic over a range of pointers.  apologies...
	  /*
	    breakpoint( const double & position, 
	    const unsigned & count,const unsigned & original_type,
	    const double & sel,const double & dom,
	    const bool & n=true) 
	  */
	  unsigned last_haplo = std::numeric_limits<unsigned>::max();
	  unsigned bpoint_found = false;
	  for( mcont_iterator itr = ibeg->mutations.begin() ;
	       itr != ibeg->mutations.end() ; ++itr )
	    {
	      last_haplo = (*itr)->hap;
	      if( (*itr)->pos <= pos )
		{
		  new_gamete1.mutations.push_back(*itr);
		}
	      else
		{
		  if(!bpoint_found)
		    {
		      list<mtype>::iterator nbitr = mutations->insert( mutations->end(),
								       breakpoint( pos,1,last_haplo,
										   0,0,true ) );
		      assert (! nbitr->checked );
		      new_gamete1.mutations.push_back(nbitr);
		      bpoint_found=true;
		    }
		  new_gamete2.mutations.push_back(*itr);
		}
	    }
	      
	  last_haplo = std::numeric_limits<unsigned>::max();
	  bpoint_found=false;
	  for( mcont_iterator itr = jbeg->mutations.begin() ; 
	       itr != jbeg->mutations.end() ; ++itr )
	    {
	      last_haplo = (*itr)->hap;
	      if( (*itr)->pos <= pos )
		{
		  new_gamete2.mutations.push_back(*itr);
		}
	      else
		{
		  if(!bpoint_found)
		    {
		      list<mtype>::iterator nbitr = mutations->insert( mutations->end(),
								       breakpoint( pos,1,last_haplo,
										   0,0,true ) );
		      assert (! nbitr->checked );
		      new_gamete2.mutations.push_back(nbitr);
		      bpoint_found=true;
		    }
		  new_gamete1.mutations.push_back(*itr);
		}
	    }
	      
	  for( mcont_iterator itr = ibeg->smutations.begin() ;
	       itr != ibeg->smutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete1.smutations.push_back(*itr);
		}
	      else
		{
		  new_gamete2.smutations.push_back(*itr);
		}
	    }
	      
	  for( mcont_iterator itr = jbeg->smutations.begin() ; 
	       itr != jbeg->smutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete2.smutations.push_back(*itr);
		}
	      else
		{
		  new_gamete1.smutations.push_back(*itr);
		}
	    }
	  std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	  std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	  std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	  std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	      
	  vtype_iterator itr = std::find(gametes->begin(),gametes->end(),new_gamete1);
	  if( itr == gametes->end() )
	    {
	      gametes->push_back(new_gamete1);
	    }
	  else
	    {
	      itr->n++;
	    }
	  itr=find(gametes->begin(),gametes->end(),new_gamete2);
	  if( itr == gametes->end() )
	    { 
	      gametes->push_back(new_gamete2);
	    }
	  else
	    {
	      itr->n++;
	    }	
	}
      NRECS--;
#ifndef NDEBUG
      NRECS_DONE++;
#endif
    }
  assert(NRECS_DONE==doublecheck);
#ifndef NDEBUG
  unsigned sum=0;
  for(ibeg=gametes->begin();ibeg!=gametes->end();++ibeg)
    {
      sum+=ibeg->n;
    }
  assert(sum==twoN);
#endif
  return doublecheck;
}

void write_output(const SimData & ohaps,
		  const vector< vector<gtype> > & vgametes,
		  const vector< list<mtype> > & vmutations,
		  const unsigned & twoNcoal,
		  const unsigned & N,
		  const unsigned & nreps_sel,
		  const unsigned & nreps_control,
		  const double & ssite_pos,
		  const char * base,
		  const unsigned & gens);

template<typename muttype>
struct sortpos2 : public binary_function<muttype,muttype,bool>
{
  inline bool operator()(const muttype & l, const muttype & r) const
  {
    return l.pos < r.pos;
  }
};

int main(int argc, char ** argv)
{
  int argument=1;
  const char * msfile = argv[argument++];
  const unsigned maxgams = atoi(argv[argument++]);
  const unsigned N = atoi(argv[argument++]);
  const double mu = 0.;
  const double s = atof(argv[argument++]);
  const double h = atof(argv[argument++]);
  const double littler = atof(argv[argument++])/2.;
  const unsigned ngens_evolve1 = atoi(argv[argument++]);
  const unsigned ngens_evolve2 = atoi(argv[argument++]);
  const unsigned ngens_evolve3 = atoi(argv[argument++]);
  const unsigned nreps_sel=atoi(argv[argument++]);
  const unsigned nreps_control=atoi(argv[argument++]);
  const char * outfile_basename = argv[argument++];
  const unsigned seed = atoi(argv[argument++]);

  vector<unsigned>ngens;
  ngens.push_back(ngens_evolve1);
  ngens.push_back(ngens_evolve2);
  ngens.push_back(ngens_evolve3);
  sort(ngens.begin(),ngens.end());

  cout << '#';
  for(int i=0;i<argc;++i)
    {
      cout << argv[i]<<' ';
    }
  cout << '\n';

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  filtering_istream in;
  in.push(gzip_decompressor());
  in.push(file_source(msfile,ios_base::in|ios_base::binary));
  if(!in)
    {
      cerr << msfile << " could not be opened for reading\n";
      exit(10);
    }

  
  vector<gtype> igametes;
  list<mtype> imutations;
  const double CRITPOS = 0.5;
  
  SimData *mstemp=new SimData(),starting_haps,ms_selsite;
  in >> *mstemp >> ws;
  in.pop();
  in.pop();
  //now, copy the desired # of sequences into starting_haps
  starting_haps.assign( &*mstemp->pbegin(),mstemp->numsites(),
			&*mstemp->begin(), maxgams );
  delete mstemp;
  RemoveInvariantColumns(&starting_haps);

  /*
    an ms block often has a problem where two adjacent sites
    are recorded as having the same position.  This is an 
    error due to not printing the output with sufficient precision
    and using ASCI files.
  */
  for(SimData::pos_iterator i = starting_haps.pbegin() ; 
      i != starting_haps.pend();++i)
    {
      if(i != starting_haps.pend()-2)
	{
	  if(*(i+1) == *(i+2))
	    {
	      //new position value is 1/2 way between i and i+1.
	      //std::cerr << *i << ' ' << *(i+1) << ' ' <<*(i+2) << "->";
	      *(i+1) = ( (*i) + ( (*(i+1)-*i)/2. ) );
	      //std::cerr << *i << ' ' << *(i+1) << ' ' << *(i+2) << "\n";
	    }
	}
    }
  //find the selected site, and initialze ms_selsite
  double mindist = std::numeric_limits<double>::max();
  bool foundmut=false;
  SimData::const_site_iterator siter,selsiteitr=starting_haps.send();  
  double ssite_pos = 1.,ssite_freq=-1;
  for( SimData::const_site_iterator itr = starting_haps.sbegin();
       !foundmut && itr != starting_haps.send(); ++itr )
    {
      double dist =  abs(itr->first - CRITPOS);
      double freq = double(count(itr->second.begin(),itr->second.end(),'1'))/double(maxgams);
      if(dist <= mindist && freq >= 0.01 && freq < 1.)
	{
	  selsiteitr = itr;
	  mindist = dist;
	  ssite_freq = freq;
	  ssite_pos = itr->first;
	}
      if(itr->first > CRITPOS && ssite_freq < 1.)
	{
	  foundmut=true;
	}
    }
  assert(ssite_freq >= 0. && ssite_freq < 1.);
  ms_selsite.assign(selsiteitr,selsiteitr+1);
  assert(ms_selsite.numsites() == 1);
  //this is a selected mutation at position ssite_pos, belonging to haplotype  numeric_limits<unsigned>::max()
  list<mtype>::iterator ssite_itr = imutations.insert(imutations.end(), breakpoint( ssite_pos, count( ms_selsite.sbegin()->second.begin(),
												      ms_selsite.sbegin()->second.end(),'1' ),
										    numeric_limits<unsigned>::max(),s,h,false ) );
  unsigned ns=0;
  for( unsigned i = 0 ; i < maxgams ; ++i)
    {
      igametes.push_back(gtype(1));
      //insert a neutral mutation that is haplotype i up to position 1.
      /*
	Note:
	We are implicitly in an infinitely-many sites world, and valid_copy works assuming
	that all positions are unique.
	Therefore, we will assign the position of each of these terminal mutations a 
	value of 1 + double(i), so that they are >= than the rightmost allowed position.
      */
      list<mtype>::iterator itr = imutations.insert(imutations.end(),breakpoint( 1.+double(i), 1,
										 i,0.,0.,true ) );
      igametes[i].mutations.push_back( itr );
      if( ms_selsite[i][0] == '1' ) //then this haplo contains selected site
	{
	  igametes[i].smutations.push_back(ssite_itr);
	  ++ns;
	}
    }

  assert( igametes.size() == maxgams );
  assert(imutations.size() == 1 + maxgams);
#ifndef NDEBUG
  for( unsigned i = 0 ; i < maxgams ; ++i)
    {
      assert( ! igametes[i].mutations.empty() );
    }
#endif
  unsigned twoNcoal=0;
  for(unsigned i=0;i<igametes.size();++i)
    {
      twoNcoal += igametes[i].n;
    }
  assert(twoNcoal == maxgams);

  //cerr << ns << ' ' << ssite_itr->n << ' ' << igametes.size() << ' ' << maxgams << ' ' << ms_selsite.size() << ' ' << twoNcoal << '\n' << ms_selsite << '\n'; exit(10);
  vector< vector<gtype> > vgametes(nreps_sel+nreps_control);
  vector< list<mtype> > vmutations(nreps_sel+nreps_control);
  vector< vector<mtype> > vfixations(nreps_sel+nreps_control);
  vector<vector<unsigned> > vftimes(nreps_sel+nreps_control);

  for( unsigned i=0;i<nreps_sel+nreps_control;++i )
    {
      valid_copy(igametes,imutations,vgametes[i],vmutations[i]);
#ifndef NDEBUG
      for( unsigned j = 0 ; j < vgametes[i].size() ; ++j)
	{
	  assert( ! vgametes[i][j].mutations.empty() );
	}
#endif
      //bottleneck the replicate
      double wbar = sample_diploid(r,&vgametes[i],twoNcoal,2*N,boost::bind(no_selection<vector<gtype>::const_iterator >,_1,_2));
      assert(check_sum(vgametes[i],2*N));
      remove_fixed_lost(&vmutations[i],&vfixations[i],&vftimes[i],1,2*N);
      //remove_fixed_lost(&vmutations[i],1,2*N);
    }

  imutations.erase(imutations.begin(),imutations.end());
  igametes.erase(igametes.begin(),igametes.end());
   vector<unsigned> gen(nreps_sel+nreps_control,0u);
  double wbar;
  for (unsigned LENGTH = 0 ; LENGTH < ngens.size() ; ++LENGTH)
    {
      if( ngens[LENGTH] )
	{
	  for( unsigned rep=0;rep<nreps_sel+nreps_control;++rep )
	    {
	      //now, do the WF sampling
	      for(;gen[rep]<ngens[LENGTH];++gen[rep])
		{
		  assert( mutations_are_sorted( vgametes[rep] ) );
		  if(rep<nreps_sel)
		    {
		      assert( nmuts_not_empty( vgametes[rep] ) );
#ifndef NDEBUG
		      for(unsigned gam = 0 ; gam < vgametes[rep].size() ; ++gam )
			{
			  for(unsigned m=0;m<vgametes[rep][gam].mutations.size();++m)
			    {
			      assert ( vgametes[rep][gam].mutations[m]->n >= vgametes[rep][gam].n );
			    }
			}
		      list<mtype> mtemp;
		      vector<gtype> gtemp;
		      valid_copy(vgametes[rep],vmutations[rep],gtemp,mtemp);
#endif
		      //wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(multiplicative_diploid(),_1,_2));
		      wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(multiplicative_diploid(),_1,_2),
					     boost::bind(update_gamete_list2(),_1,_2));
		      
#ifndef NDEBUG
		      if (! nmuts_not_empty(vgametes[rep]) )
			{
			  for(unsigned gam=0 ; gam<vgametes[rep].size() ; ++gam )
			    {
			      if( vgametes[rep][gam].mutations.empty() )
				{
				  cout << gtemp[gam].mutations.size() << '\t' << gtemp[gam].smutations.size() << ' '
				       << vgametes[rep][gam].mutations.size() << '\t' << vgametes[rep][gam].smutations.size() << '\n';
				  for( unsigned i = 0 ; i < gtemp[gam].mutations.size() ; ++i )
				    {
				      cout << gtemp[gam].mutations[i]->n << ' ';
				    }
				  cout << "|";
				  for( unsigned i = 0 ; i < gtemp[gam].smutations.size() ; ++i )
				    {
				      cout << gtemp[gam].smutations[i]->n << ' ';
				    }
				  cout << '\n';
				  for(unsigned i=0;i<vgametes[rep][gam].mutations.size();++i)
				    {
				      cout << vgametes[rep][gam].mutations[i]->n << '\n';
				    }
				  cout << "|";
				  for(unsigned i=0;i<vgametes[rep][gam].smutations.size();++i)
				    {
				      cout << vgametes[rep][gam].smutations[i]->n << ' ';
				    }
				  cout << '\n';
				}
			    }
			}
#endif
		      assert( nmuts_not_empty( vgametes[rep] ) );
		    }
		  else
		    {
		      assert( nmuts_not_empty( vgametes[rep] ) );
		      //wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(no_selection<vector<gtype>::const_iterator>,_1,_2));
		      wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(no_selection<vector<gtype>::const_iterator>,_1,_2),
					     boost::bind(update_gamete_list2(),_1,_2));
		      assert( nmuts_not_empty( vgametes[rep] ) );
		    }
		  assert( nmuts_not_empty( vgametes[rep] ) );
		  //remove_fixed_lost(&vmutations[rep],&vfixations[rep],&vftimes[rep],gen[rep],2*N);
		  remove_lost(&vmutations[rep],gen[rep],2*N);
		  assert( nmuts_not_empty( vgametes[rep] ) );
		  assert(check_sum(vgametes[rep],2*N));
		  unsigned nrec = recombine_breakpoints(r, &vgametes[rep], &vmutations[rep],2*N, littler, boost::bind(gsl_rng_uniform,r));
		  assert(check_sum( vgametes[rep],2*N ) );
		  assert( mutations_are_sorted( vgametes[rep] ) );
		  assert( nmuts_not_empty( vgametes[rep] ) );
		  remove_lost(&vmutations[rep],gen[rep],2*N);
		  //remove_fixed_lost(&vmutations[rep],&vfixations[rep],&vftimes[rep],gen[rep],2*N);
		}
	    }
	  write_output(starting_haps, vgametes, vmutations, twoNcoal, N, nreps_sel, nreps_control, ssite_pos, outfile_basename, gen[0]);
	}
    }
}

void write_output(const SimData & ohaps,
		  const vector< vector<gtype> > & vgametes,
		  const vector< list<mtype> > & vmutations,
		  const unsigned & twoNcoal,
		  const unsigned & N,
		  const unsigned & nreps_sel,
		  const unsigned & nreps_control,
		  const double & ssite_pos,
		  const char * base,
		  const unsigned & gens)
/*
  Take 2!
*/
{
   ostringstream o;
   o << base << '.' << gens << ".gz";
   filtering_ostream out;
   out.push(gzip_compressor());
   out.push(file_sink(o.str().c_str(),ios_base::binary|ios_base::out));

   o.str(string());
   double ddummy;
   for(SimData::const_site_iterator i = ohaps.sbegin() ; i != ohaps.send() ; ++i)
     {
       //get indexes of which founder haps have this mut
       string::const_iterator sci = find(i->second.begin(),i->second.end(),'1');
       vector<unsigned> hindexes;
       while(sci != i->second.end())
	 {
	   hindexes.push_back( sci - i->second.begin());
	   sci = find(sci + 1,i->second.end(),'1');
	 }
       vector< double > repfreqs(nreps_sel + nreps_control,0.);
       for(unsigned rep=0;rep < (nreps_sel + nreps_control) ; ++rep)
	 {
	   for(unsigned gam = 0 ; gam < vgametes[rep].size() ; ++gam)
	     {
	       bool found = false;
	       for(unsigned m=0;!found&&m<vgametes[rep][gam].mutations.size();++m)
		 {
		   double lastpos = (m > 0) ? vgametes[rep][gam].mutations[m-1]->pos : 0.;
		   if(vgametes[rep][gam].mutations[m]->pos >= i->first && i->first > lastpos)
	       		 {
	       		   if(find(hindexes.begin(),hindexes.end(),vgametes[rep][gam].mutations[m]->hap) != hindexes.end() )
	       		     {
	       		       found = true;
	       		       repfreqs[rep] +=  vgametes[rep][gam].n;
	       		     }
	       		 }
		 }
	     }
	 }
       transform(repfreqs.begin(),repfreqs.end(),repfreqs.begin(),bind2nd(divides<double>(),double(2*N)));
       o << i->first << '\t'
	   << double(count(i->second.begin(),i->second.end(),'1'))/double(twoNcoal) << '\t'
	   << (i->first == ssite_pos) << '\t';
       for(unsigned j=0;j<repfreqs.size();++j)
	 {
	   o << repfreqs[j];
	   if(j<repfreqs.size()-1)
	     {
	       o << '\t';
	     }
	   else
	     {
	       o << '\n';
	     }
	 }
       if( o.str().size() > 10000000 )
	 {
	   out << o.str();
	   o.str(string());
	 }
     }
   if(o.str().size() > 0)
     {
       out << o.str();
     }
   out.pop();
   out.pop();
}
