#ifndef __FWDPP_SUGAR_MUTATION_INFSITES_HPP__
#define __FWDPP_SUGAR_MUTATION_INFSITES_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/sugar/popgenmut.hpp>

namespace KTfwd
{
  /*!
    \brief Generic function object implementing 
    the infinitely-many sites mutation model for 
    "standard" population-genetic scenarios
    \ingroup sugar
   */
  struct infsites
  {
    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename queue_t,
	     typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,
				   typename mlist_t::iterator>::type
    operator()(queue_t & recycling_bin,
	       mlist_t * mutations,
	       gsl_rng * r, lookup_table_t * lookup,
	       const uint_t & generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      static_assert(std::is_same<typename mlist_t::value_type,KTfwd::popgenmut>::value,
		    "mlist_t::value_type must be KTfwd::mutation");
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      bool selected = (gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate));
      if (!recycling_bin.empty())
	{
	  auto i = recycling_bin.front();
	  i->pos=pos;
	  i->s=(selected) ? smaker() : 0.;
	  i->h=(selected) ? hmaker() : 0.;
	  i->g=generation;
	  i->n=1u;
	  recycling_bin.pop();
	  return i;
	}
      return mutations->emplace(mutations->end(),pos,(selected)?smaker():0.,(selected)?hmaker():0.,generation,1u);
    }

    /*!
      \brief Overload for different position distributions for neutral and non-neutral variants

      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param nposmaker A policy that returns the position of new neutral mutation
      \param sposmaker A policy that returns the position of new non-neutral mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename queue_t,
	     typename mlist_t,
	     typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,
				   typename mlist_t::iterator>::type
    operator()(queue_t & recycling_bin,
	       mlist_t * mutations,
	       gsl_rng * r, lookup_table_t * lookup,
	       const uint_t & generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  if(!recycling_bin.empty())
	    {
	      auto rv = recycling_bin.front();
	      recycling_bin.pop();
	      rv->pos=pos;
	      rv->s = smaker();
	      rv->h = hmaker();
	      rv->n = 1u;
	      rv->g = generation;
	      return rv;
	    }
	  return mutations->emplace(mutations->end(),pos,smaker(),hmaker(),generation,1);
	}
      //Establish position of new mutation
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      if(!recycling_bin.empty())
	{
	  auto rv = recycling_bin.front();
	  recycling_bin.pop();
	  rv->pos=pos;
	  rv->s=0.;
	  rv->h=0.;
	  rv->g=generation;
	  rv->n=1;
	  return rv;
	}
      return mutations->emplace(mutations->end(),pos,0.,0.,generation,1);
    }

    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename queue_t,
	     typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,
				   typename mlist_t::iterator>::type
    operator()(queue_t & recycling_bin,
	       mlist_t * mutations,
	       gsl_rng * r, lookup_table_t * lookup,
	       const uint_t * generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      static_assert(std::is_same<typename mlist_t::value_type,KTfwd::popgenmut>::value,
		    "mlist_t::value_type must be KTfwd::mutation");
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      bool selected = (gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate));
      if (!recycling_bin.empty())
	{
	  auto i = recycling_bin.front();
	  i->pos=pos;
	  i->s=(selected) ? smaker() : 0.;
	  i->h=(selected) ? hmaker() : 0.;
	  i->g=*generation;
	  i->n=1u;
	  recycling_bin.pop();
	  return i;
	}
      return mutations->emplace(mutations->end(),pos,(selected)?smaker():0.,(selected)?hmaker():0.,*generation,1u);
    }

    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename queue_t,
	     typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,mutation>::value,
				   typename mlist_t::iterator>::type
    operator()(queue_t & mutation_recycling_bin,
	       mlist_t * mutations,
	       gsl_rng * r, lookup_table_t * lookup,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      static_assert(std::is_same<typename mlist_t::value_type,KTfwd::mutation>::value,
		    "mlist_t::value_type must be KTfwd::mutation");
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      bool selected = (gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate));
      if (!mutation_recycling_bin.empty())
	{
	  auto i = mutation_recycling_bin.front();
	  i->pos=pos;
	  i->s=(selected) ? smaker() : 0.;
	  i->h=(selected) ? hmaker() : 0.;
	  i->n=1u;
	  mutation_recycling_bin.pop();
	  return i;
	}
      return mutations->emplace(mutations->end(),pos,(selected)?smaker():0.,(selected)?hmaker():0.,1u);
    }

    /*!
      \brief Overload for different position distributions for neutral and non-neutral variants

      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param nposmaker A policy that returns the position of new neutral mutation
      \param sposmaker A policy that returns the position of new non-neutral mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
    */
    template<typename queue_t,
	     typename mlist_t,
	     typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,mutation>::value,
				   typename mlist_t::iterator>::type
    operator()(queue_t & recycling_bin,
	       mlist_t * mutations,
	       gsl_rng * r, lookup_table_t * lookup,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if ( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  if(!recycling_bin.empty())
	    {
	      auto rv = recycling_bin.front();
	      recycling_bin.pop();
	      rv->pos=pos;
	      rv->s=smaker();
	      rv->h=hmaker();
	      rv->n=1u;
	      return rv;
	    }
	  return mutations->emplace(mutations->end(),pos,smaker(),1,hmaker());
	}
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      //return a neutral mutation
      if(!recycling_bin.empty())
	{
	  auto rv = recycling_bin.front();
	  recycling_bin.pop();
	  rv->pos=pos;
	  rv->s=0.;
	  rv->h=0.;
	  rv->n=1u;
	  return rv;
	}
      return mutations->emplace(mutations->end(),pos,0.,1,0.);
    }
  };
}

#endif
