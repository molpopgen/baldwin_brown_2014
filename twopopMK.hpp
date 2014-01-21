#ifndef __TWOPOPMK_HPP__
#define __TWOPOPMK_HPP__

#include <diploid.hh>
#include <Sequence/SimData.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/PolySIM.hpp>
#include <ctime>
#include <iostream>
#include <fstream>

struct mutation_with_age : public mutation_base
{
  mutable unsigned o;
  double s,h;
  char label;
  mutation_with_age( const double & position, const double & sel_coeff,
		     const unsigned & count,const unsigned & origin, const char & ch,
		     const double & dominance = 0.5, const bool & n=true); 
  bool operator==(const mutation_with_age & rhs) const;
};

typedef std::list<mutation_with_age > mlist;
typedef gamete_base<mutation_with_age, mlist> gtype;
typedef std::vector<gtype> gvector;
typedef std::vector<mutation_with_age> mvector;

struct samples
{
  std::vector< std::pair<double,std::string> > A,S;
  mvector Amuts,Smuts;
  unsigned FR,FS,PR,PS;
  samples( const  std::vector< std::pair<double,std::string> > & a, const std::vector< std::pair<double,std::string> > & s,
	   const mvector & amuts, const mvector & smuts,
	   const unsigned & fr, const unsigned &fs,
	   const unsigned & pr, const unsigned &ps ) :
    A(a),S(s),Amuts(amuts),Smuts(smuts),
    FR(fr),FS(fs),PR(pr),PS(ps)
  {
  }
};

struct sortpos2 : public std::binary_function<mutation_with_age,mutation_with_age,bool>
{
  inline bool operator()(const mutation_with_age & a, const mutation_with_age & b) const
  {
    return (a.pos < b.pos);
  }
};

inline double simple_map(gsl_rng * r, const bool & linked)
{
  if (linked)
    return (gsl_ran_flat(r,-10.,10.));
  
  return gsl_ran_flat(r,0.,1.);
}

samples ms_sample(gsl_rng * r,
		  const gvector & gametes1,
		  const gvector & gametes2,
		  const mlist & mdaughter1,
		  const mlist & mdaughter2,
		  const mvector & fixations1,
		  const mvector & fixations2,
		  const unsigned & n1, const unsigned & n2,const unsigned & N1, const unsigned & N2,
		  const bool & linked);

void write_output(gsl_rng *r,  const gvector & gdaughter1,const  gvector & gdaughter2,
		  const mlist & mdaughter1, const mlist & mdaughter2,
		  const mvector & fdaughter1, const mvector & fdaughter2,
		  const std::vector<unsigned> & ftimes1, const std::vector<unsigned> & ftimes2,
		  const unsigned & samplesize, const unsigned & samplesize2, const unsigned &N0,
		  const char * msfile, const char * mutfile,const char * mkfile,const bool & linked);

#endif
