#include <diploid.hh>
#include <initms.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <isbinary.hpp>
#include <sstream>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <Sequence/Recombination.hpp>

using namespace std;
using namespace Sequence;
using namespace boost::iostreams;

#ifndef NDEBUG
bool check_sum(const std::vector<gamete> & gametes, const unsigned & twoN)
{
  unsigned check=0;
  for(std::vector<gamete>::const_iterator i=gametes.begin();i<gametes.end();++i)
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


void write_output(const map<double,unsigned> & ifreqs,
		  const vector< vector<gamete> > & vgametes,
		  const vector< list<mutation> > & vmutations,
		  const vector< vector<mutation> > & vfixations,
		  const unsigned & twoNcoal,
		  const unsigned & N,
		  const unsigned & nreps_sel,
		  const unsigned & nreps_control,
		  const double & ssite_pos,
		  const char * base,
		  const unsigned & gens);

//should be much faster
void write_output( map<double,unsigned> & ifreqs,
		   vector< map<double, list<mutation>::iterator> > & mpointers,
		   const unsigned & twoNcoal,
		   const unsigned & N,
		   const double & ssite_pos,
		   const char * base,
		  const unsigned & gens);

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

  vector<gamete> igametes;
  list<mutation> imutations;
  cerr << "initializing from coalescent sample...";
  init_with_ms_from_stream(r,in,&igametes,&imutations,maxgams);
  cerr << "done" << endl;
  unsigned twoNcoal=0;
  for(unsigned i=0;i<igametes.size();++i)
    {
      twoNcoal += igametes[i].n;
    }
  map<double,unsigned> ifreqs;
  for(list<mutation>::const_iterator i = imutations.begin();
      i != imutations.end();++i)
    {
      ifreqs[i->pos]=i->n;
    }
  const double CRITPOS = 0.5;
  //find the selected site
  double mindist = std::numeric_limits<double>::max();
  bool foundmut=false;
  double ssite_pos = 1.,ssite_freq=-1;
  for(list<mutation>::iterator itr=imutations.begin();!foundmut&&itr!=imutations.end();++itr)
    {
      double dist =  abs(itr->pos - CRITPOS);
      double freq = double(itr->n)/double(twoNcoal);
      if(dist <= mindist && freq >= 0.01 && freq < 1.)
	{
	  ssite_pos = itr->pos;
	  ssite_freq = itr->n/double(twoNcoal);
	  mindist = dist;
	}
      if( itr-> pos > CRITPOS && ssite_freq < 1. )
	{
	  foundmut = true;
	}
    }
  //make the mutation selected!
  foundmut = false;
  list<mutation>::iterator ssite_itr;
  for( list<mutation>::iterator itr=imutations.begin();!foundmut&&itr!=imutations.end();++itr)
    {
      if(itr->pos == ssite_pos )
	{
	  foundmut = true;
	  itr->s = s;
	  itr->h = h;
	  itr->neutral = false;
	  ssite_itr = itr;
	}
    }

  //update the gametes so that this one is in their smutations
  for(unsigned i=0;i<igametes.size();++i)
    {
      for( unsigned j=0;j<igametes[i].mutations.size();++j)
	{
	  if(igametes[i].mutations[j] == ssite_itr)
	    {
	      igametes[i].smutations.push_back(igametes[i].mutations[j]);
	      igametes[i].mutations.erase(igametes[i].mutations.begin()+j);
	      j = igametes[i].mutations.size();
	    }
	}
    }
  vector< vector<gamete> > vgametes(nreps_sel+nreps_control);
  vector< list<mutation> > vmutations(nreps_sel+nreps_control);
  vector< vector<mutation> > vfixations(nreps_sel+nreps_control);
  vector<vector<unsigned> > vftimes(nreps_sel+nreps_control);
  cerr << "initalizing replicate populations...";

  for( unsigned i=0;i<nreps_sel+nreps_control;++i )
    {
      valid_copy(igametes,imutations,vgametes[i],vmutations[i]);
      //bottleneck the replicate
      double wbar = sample_diploid(r,&vgametes[i],twoNcoal,2*N,boost::bind(no_selection<vector<gamete>::const_iterator >,_1,_2));
      assert(check_sum(vgametes[i],2*N));
      //remove_fixed_lost(&vmutations[i],&vfixations[i],&vftimes[i],1,2*N);
      for(unsigned j=0;j<vmutations.size();++j)
	{
	  for( list<mutation>::iterator k = vmutations[j].begin() ; 
	       k != vmutations[j].end() ; ++k )
	    {
	      k->checked = false;
	    }
	}
    }

  /*
    By using the data structure below, we can no longer call "remove_fixed_lost".  We have
    to take that performance hit in order to speed up the output, which is quite slow otherwise.
  */
  vector< map<double, list<mutation>::iterator> >  mpointers(vmutations.size());  //pointers to each mutation in each replicate
  for( unsigned rep=0;rep<vmutations.size();++rep )
    {
      for( list<mutation>::iterator i = vmutations[rep].begin() ; i != vmutations[rep].end() ; ++i )
	{
	  mpointers[rep].insert( make_pair(i->pos,i) );
	}
    }
  imutations.erase(imutations.begin(),imutations.end());
  igametes.erase(igametes.begin(),igametes.end());
  //imutations.clear();
  //igametes.clear();
  cerr << "done" <<endl;
  //now, evolve every replicate
  //vector<vector<pair<double,unsigned> > > pos(nreps_sel+nreps_control);
  //unsigned gen = 0;
  vector<unsigned> gen(nreps_sel+nreps_control,0u);
  double wbar;
  for (unsigned LENGTH = 0 ; LENGTH < ngens.size() ; ++LENGTH)
    {
      if( ngens[LENGTH] )
	{
	  cerr << "evolving until " << ngens[LENGTH] << " generations...";
	  for( unsigned rep=0;rep<nreps_sel+nreps_control;++rep )
	    {
	      //now, do the WF sampling
	      for(;gen[rep]<ngens[LENGTH];++gen[rep])
		{
		  if(rep<nreps_sel)
		    {
		      wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(multiplicative_diploid(),_1,_2));
		    }
		  else
		    {
		      wbar = sample_diploid(r,&vgametes[rep],2*N,boost::bind(no_selection<vector<gamete>::const_iterator>,_1,_2));
		    }
		  //remove_fixed_lost(&vmutations[rep],&vfixations[rep],&vftimes[rep],gen[rep],2*N);
		  for(unsigned j=0;j<vmutations.size();++j)
		    {
		      for( list<mutation>::iterator k = vmutations[j].begin() ; 
			   k != vmutations[j].end() ; ++k )
			{
			  k->checked = false;
			}
		    }
		  assert(check_sum(vgametes[rep],2*N));
		  unsigned nrec = recombine(r, &vgametes[rep], 2*N, littler, boost::bind(gsl_rng_uniform,r));
		  //assert(check_sum(gametes[rep],2*N));
		}
	    }
	  cerr << "done\noutputting...";
      
	  /*
	  write_output(ifreqs,vgametes,vmutations,vfixations,
		       twoNcoal,N,nreps_sel,nreps_control,
		       ssite_pos,outfile_basename,gen[0]);
	  */
	  write_output(ifreqs,mpointers,twoNcoal,N,ssite_pos,outfile_basename,gen[0]);
	  cerr << "done\n";
	}
    }



  cerr << "done" << endl;
}

void write_output(const map<double,unsigned> & ifreqs,
		  const vector< vector<gamete> > & vgametes,
		  const vector< list<mutation> > & vmutations,
		  const vector< vector<mutation> > & vfixations,
		  const unsigned & twoNcoal,
		  const unsigned & N,
		  const unsigned & nreps_sel,
		  const unsigned & nreps_control,
		  const double & ssite_pos,
		  const char * base,
		  const unsigned & gens)
{
 
   ostringstream o;
   o << base << '.' << gens << ".gz";
   filtering_ostream out;
   out.push(gzip_compressor());
   out.push(file_sink(o.str().c_str(),ios_base::binary|ios_base::out));
   list<mutation>::const_iterator lmi;
   vector<mutation>::const_iterator vmi;
  //  for( unsigned i = 0 ; i < markers.size() ; ++i )
  //   {
      for( map<double,unsigned>::const_iterator i = ifreqs.begin() ; 
      	   i != ifreqs.end() ; ++i )
      	{
      	  out << i->first << '\t' << double(i->second)/double(twoNcoal) << '\t'
      	      << (i->first == ssite_pos) << '\t';
      	  for(unsigned rep = 0 ; rep < nreps_sel+nreps_control ; ++rep)
      	    {
      	      lmi = find_if(vmutations[rep].begin(),
      			    vmutations[rep].end(),
      			    boost::bind(mutation_at_pos(),_1,i->first));
      	      if(lmi != vmutations[rep].end())
      		{
      		  out << double(lmi->n)/double(2*N);
      		}
      	      else
      		{
      		  vmi = find_if(vfixations[rep].begin(),
      				vfixations[rep].end(),
      				boost::bind(mutation_at_pos(),_1,i->first));
      		  if(vmi != vfixations[rep].end())
      		    {
      		      out << double(vmi->n)/double(2*N);
      		    }
      		  else
      		    {
      		      out << 0;
      		    }
      		}
      	      if(rep < nreps_sel+nreps_control-1)
      		{
      		  out << '\t';
      		}
      	    }
      	  out << '\n';
      	}
  out.pop();
  out.pop();
}

void write_output( map<double,unsigned> & ifreqs,
		   vector< map<double, list<mutation>::iterator> > & mpointers,
		   const unsigned & twoNcoal,
		   const unsigned & N,
		   const double & ssite_pos,
		   const char * base,
		   const unsigned & gens)
{
 
   ostringstream o;
   o << base << '.' << gens << ".gz";
   filtering_ostream out;
   out.push(gzip_compressor());
   out.push(file_sink(o.str().c_str(),ios_base::binary|ios_base::out));

   for( map<double,unsigned>::iterator i = ifreqs.begin() ; 
	i != ifreqs.end() ; ++i )
     {
       out << i->first << '\t' << double(i->second)/double(twoNcoal) << '\t'
	   << (i->first == ssite_pos) << '\t';
       for( unsigned rep = 0 ; rep < mpointers.size() ; ++rep )
	 {
	   out << double(mpointers[rep][i->first]->n)/double(2*N);
	   if(rep < mpointers.size()-1)
	     {
	       out << '\t';
	     }
	 }
       out << '\n';
     }
   out.pop();
   out.pop();
}
