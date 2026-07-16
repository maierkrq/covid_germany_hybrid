/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MULTIVARIATIONALUTIL_HH
#define MULTIVARIATIONALUTIL_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/include/fold.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/mpl/int.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>


#include "fem/fixdune.hh"
#include "fem/functionspace.hh"
#include "fem/linearspace.hh"
#include "io/iobase.hh" // for paddedString
#include "linalg/crsutil.hh"

/**
 * @file
 * @brief  Variables and their descriptions
 * @author Martin Weiser
 *
 * This file contains a number of core classes of Kaskade, needed to define a problem. The important ones are:
 *
 * - VariableDescription: elementary information about a single variable
 * - VariableSetDescription: information about a set of variables and associated FEFunctionSpace s
 * - VariableSet: a set of variables with all the data stored in FunctionSpaceElement s
 *
 * The other classes are auxilliary
 */


//---------------------------------------------------------------------

namespace Kaskade
{
  /**
   * \brief A class storing elementary information about a single variable.
   * \ingroup variables
   * This includes the number of components, the id and name of the
   * variable, and the index of the associated FEFunctionSpace from
   * which the variable comes.
   *
   * \tparam spaceId the index of the FEFunctionSpace to which the variable belongs.
   *   This references an external container of FE spaces.
   * \tparam components the number of components that this variable has
   *   (1 indicates a scalar variable, values greater than 1 a vectorial
   *   variable). Has to be greater than 0.
   * \tparam Id the index of the variable by which it can be retrieved from a
   *   boost::fusion container of variables.
   *
   * Note that the order of template parameters is important. A more convenient interface
   * is provided by \ref Variable, where the order of template parameters can be arbitrary.
   */
  template <int spaceId, int components, int Id=-1>
  struct VariableDescription
  {
    /// number of this variable
    static int const id = Id;
    /// number of component of this variable
    static int const m = components;
    /// number of the space, this variable is associated with
    static int const spaceIndex = spaceId;
  };


  /**
   * \brief Helper class for specifying the FE space index of a variable.
   * 
   * This is intended to be used as template argument for \ref Variable.
   */
  template <int spaceIndex>
  struct SpaceIndex : public VariableDescription<spaceIndex,-1,-1> {};

  /**
   * \brief Helper class for specifying the number of components of a variable.
   * 
   * This is intended to be used as template argument for \ref Variable.
   */
  template <int components>
  struct Components : public VariableDescription<-1,components,-1> {};

  /**
   * \brief DEPRECATED
   */
  template <int id>
  struct [[deprecated]] VariableId : public VariableDescription<-1,-1,id> {};


  /**
   * \brief A class defining elementary information about a single variable.
   * \ingroup variables
   *
   * This includes the number of components, the name of the
   * variable, and the index of the associated FEFunctionSpace from
   * which the variable comes.
   *
   * This class differs from \ref VariableDescription mainly in how it is created. The template
   * parameters should be \ref SpaceIndex and \ref Components, in arbitrary order (but one of each).
   *
   * Example (e.g. for a Stokes setting):
   * \code
   * using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<3>>,
   *                                                    Variable<Components<1>,SpaceIndex<1>>>;
   * \endcode
   */
  template <typename A, typename B, typename C = VariableDescription<-1, -1, -1>>
  struct Variable 
  {
    static int const spaceIndex = A::spaceIndex>=0? A::spaceIndex : B::spaceIndex>=0? B::spaceIndex : C::spaceIndex;
    static int const m  = A::m>=0? A::m : B::m>=0? B::m : C::m;

    static_assert(spaceIndex>=0,"Space index has to be nonnegative (space indices cover a contiguous range from 0 to k-1).");
    static_assert(m>=0,"Number of components has to be nonnegative (in fact, 0 makes rarely sense as well).");

    /**
     * \brief Constructor.
     * \param name The name of the variable. This is used for output.
     */
    Variable(std::string const& name_)
    : name(name_) {}
    
    /**
     * \brief Constructor.
     * This will create a default name based on the variable id.
     */
    Variable() = default;
    
    std::string name;
  };
  

  //---------------------------------------------------------------------

  /**
   * \brief A boost::fusion functor for generating function space elements for given variables.
   *
   * \tparam Spaces: a boost::fusion container of pointers to FEFunctionSpace s.
   */
  template <class Spaces>
  struct ConstructElement
  {
    ConstructElement(Spaces const& spaces_):
      spaces(spaces_)
    {}

    /**
     * \brief For a given variable, defines the correct FunctionSpaceElement type. 
     * 
     * This depends on the spaceIndex and the number of components defined by the Variable (which is a VariableDescription).
     */
    template <typename Arg> struct result {};

    template <class Variable>
    struct result<ConstructElement<Spaces>(Variable)>
    {
      typedef typename boost::remove_reference<Variable>::type Var;
      typedef typename SpaceType<Spaces,Var::spaceIndex>::type Space;
      typedef typename Space::template Element<Var::m>::type type;
    };

    /**
     * \brief Returns a FunctionSpaceElement element from the FEFunctionSpace associated to the Variable.
     */
    template <class Variable>
    typename result<ConstructElement<Spaces>(Variable)>::type operator()(Variable const& variable) const
    {
      return typename result<ConstructElement<Spaces>(Variable)>::type(*boost::fusion::at_c<Variable::spaceIndex>(spaces));
    }

  private:
    Spaces const& spaces;
  };

  /**
   * \brief A boost::fusion functor for generating coefficient vectors for given variables.
   *
   * \tparam Spaces: a boost::fusion container of pointers to FEFunctionSpace s.
   */
  template <class Spaces>
  struct ConstructCoefficientVector
  {
    ConstructCoefficientVector(Spaces const& spaces_):
      spaces(spaces_)
    {}

    /**
     * For a given variable (VariableDescription), defines the correct
     * coefficient vector type. This depends on the space index and the
     * number of components defined by the variable.
     */
    template <typename Arg> struct result {};

    template <class Variable>
    struct result<ConstructCoefficientVector<Spaces>(Variable)>
    {
      typedef typename boost::remove_reference<Variable>::type Var;
      typedef typename SpaceType<Spaces,Var::spaceIndex>::type Space;
      typedef typename Space::template Element<Var::m>::type::StorageType StorageType;

      typedef Dune::BlockVector<typename StorageType::block_type> type;
    };

    /**
     * \brief Returns a FunctionSpaceElement from the FEFunctionSpace associated to the variable.
     *
     * This is currently not very efficient. This may change when Dune starts supporting
     * move-semantics(C++11). Alternatively it may make sense to return only the dof instead
     * of the full matrix and construct the latter later on.
     */
    template <class Variable>
    typename result<ConstructCoefficientVector<Spaces>(Variable)>::type operator()(Variable const& variable) const
    {
      return typename result<ConstructCoefficientVector<Spaces>(Variable)>::type(boost::fusion::at_c<Variable::spaceIndex>(spaces)->degreesOfFreedom());
    }

  private:
    Spaces const& spaces;
  };

  //---------------------------------------------------------------------
  
  /**
   * \cond internals
   */
  namespace Variables_Detail 
  {
    using namespace boost::fusion;
    
    // A functor constructing representations of consecutive subsets of variables
    template <class Variables, class RepresentationConstructor,
              int first=0, int last=result_of::size<Variables>::type::value>
    class VariableRangeCreator
    {
      typedef typename result_of::begin<Variables>::type       Begin;
      typedef typename result_of::advance_c<Begin,first>::type First;
      typedef typename result_of::advance_c<Begin,last>::type  Last;
      
    public:
      typedef typename boost::fusion::transform_view<iterator_range<First,Last> const,
                                                     RepresentationConstructor> View;
      using type = typename result_of::as_vector<View>::type;
      
      static View apply(RepresentationConstructor const& c)
      {
        using namespace boost::fusion;
        Variables vars;
        return transform(iterator_range<First,Last>(advance_c<first>(begin(vars)),advance_c<last>(begin(vars))),c);
      }
    };

    // Metafunction class for boost mpl transform. It takes a variable description and an id,
    // and overwrites the variable id by the given id.
    struct AddIdToVariableDescription 
    {
      template <class Vd, class Int>
      struct apply 
      {
        using type = Kaskade::VariableDescription<Vd::spaceIndex, Vd::m, Int::value>;
      };
    };

    // Given a type of a sequence of VariableDescriptions a new sequence type of VariableDescriptions is created 
    // where the id of each VariableDescription is the position in the sequence
    template <typename IncompleteVariableDescriptions>
    struct ImplicitIdVariableDescriptions 
    {
      using size = typename result_of::size<IncompleteVariableDescriptions>::type;
      using type = typename result_of::as_vector<typename boost::mpl::transform<IncompleteVariableDescriptions, 
                                                                                boost::mpl::range_c<int,0,size::value>,
                                                                                AddIdToVariableDescription>::type
                                                 >::type;
    };
    
    // Given a sequence of variables and a sequence of spaces, returns the minimal 
    // continuity of all variables.
    template <typename Variables, typename Spaces>
    class Continuity
    {
      struct Max 
      {
        template <typename C, typename Var>
        auto operator()(C,Var) const 
        {
          using SpacePtr = typename result_of::value_at_c<Spaces,Var::spaceIndex>::type;
          constexpr int m = std::remove_pointer_t<SpacePtr>::continuity;
          return typename boost::mpl::min<C,boost::mpl::int_<m>>::type();
        }
      };
      
    public:
      static int const value = result_of::fold<Variables,boost::mpl::int_<1000>,Max>::type::value;
    };
    
  }
  /**
   * \endcond
   */

  //---------------------------------------------------------------------

  /**
   * \ingroup variables
   * \brief A class for storing a heterogeneous collection of FunctionSpaceElement s.
   * \tparam Descriptions a VariableSetDescription class 
   *
   * Basic vector operations and I/O are supported.
   *
   * The type of a VariableSet is defined via a VariableSetDescription, which contains information
   * about the variables in form of a boost::fusion::vector of FEFunctionSpace s and VariableDescription s.
   * The type of a VariableSet can be determined by VariableSetDescription::VariableSet.
   *
   * VariableSet s contains a collection of FunctionSpaceElement s (vars) and a reference to its VariableSetDescription (descriptions). 
   * Simultaneous evaluation of all functions is supported by \ref evaluateVariables.
   *
   * Usually, solutions are stored in VariableSets. Access to data is via the public member
   * variable data, which is a boost::fusion::vector of \ref FunctionSpaceElement s
   */
  template <class VSDescriptions>
  class VariableSet: public LinearProductSpace<typename VSDescriptions::Scalar, 
                                               typename VSDescriptions::RepresentationData>
  {
  public:
    /// Type of the VariableSetDescription
    typedef VSDescriptions Descriptions;

    /// boost::fusion::vector of data elements (of type FunctionSpaceElement)
    typedef typename VSDescriptions::RepresentationData Functions;

  private:
    typedef VariableSet<Descriptions> Self;
    typedef typename SpaceType<typename Descriptions::Spaces,0>::type::Scalar Scalar;
    typedef LinearProductSpace<Scalar,Functions> Base;

    /// Grid type
    typedef typename Descriptions::Grid Grid;
    using GridView = typename Descriptions::GridView;

  public:
    typedef LinearProductSpace<Scalar,Functions> LinearSpace;

    /// Copy constructor
    VariableSet(Self const& vs):
      Base(vs.data), descriptions(vs.descriptions)
    {}

    /**
     * \brief Constructor
     * 
     * This creates a zero-initialized set of variables.
     * \param d the VariableSetDescription object (has to exist during the lifetime of this VariableSet).
     */
    explicit VariableSet(Descriptions const& d):
        Base(Variables_Detail::VariableRangeCreator<typename Descriptions::Variables,ConstructElement<typename Descriptions::Spaces> >::apply(
          ConstructElement<typename Descriptions::Spaces>(d.spaces))), 
        descriptions(d)
    {}

    /// Assignment
    // TODO: The reference to the variable set description is not modified. Danger is,
    // that the assigned-to function space elements are now elements from different spaces
    // than before (and stated by the variable set description). Maybe switch to a (private)
    // pointer member and provide a constant getter?
    Self& operator=(Self const& v)
    {
      if (this!=&v)
        this->data = v.data;
      return *this;
    }

    /**
     * \brief Assignment from other (compatible) types.
     *
     * Compatible types are those for which the base class LinearProductSpace<...>
     * has a matching assignment operator. This includes, in particular, the coefficient
     * vector representation Descriptions::CoefficientVectorRepresentation::type, i.e. a
     * LinearProductSpace<Scalar,DuneBlockVector<...>> with matching entries and lengths.
     */
    template <class S>
    Self& operator=(S const& s)
    {
      static_cast<Base&>(*this) = s;
      return *this;
    }
    
    /**
     * \name Evaluation of variables
     * @{
     */
    
    /**
     * \brief Value of an individual variable in the set.
     * \tparam i the index of the variable to evaluate 
     * \tparam Evaluators a boost::fusion sequence of evaluators, corresponding to 
     *         the list of spaces of this variable set
     */
    template <int i, class Evaluators>
    auto value(Evaluators const& evals) const 
    {
      static_assert(boost::fusion::result_of::size<Evaluators>::type::value==boost::fusion::result_of::size<typename Descriptions::Spaces>::type::value,
                    "Evaluators have to match the space list");
      constexpr int spaceIndex = Descriptions::template spaceIndex<i>;
      return component<i>(*this).value(boost::fusion::at_c<spaceIndex>(evals));
    }

    /**
     * \brief Value of an individual variable in the set.
     * \tparam i the index of the variable to evaluate
     * \param cell the cell in which to evaluate the variable
     * \param xi the local coordinate within the cell at which to evaluate
     */
    template <int i>
    auto value(Cell<GridView> const& cell, LocalPosition<GridView> const& xi) const
    {
      return component<i>(*this).value(cell,xi);
    }
    
    /**
     * \brief Value of an individual variable in the set. 
     * \tparam i the index of the variable to evaluate
     * \param x the global coordinate within the geometry at which to evaluate
     * \WARNING This method is comparatively DEAD SLOW.  
     */
    template <int i>
    auto value(GlobalPosition<GridView> const& x) const
    {
      return component<i>(*this).value(x);
    }

    /**
     * \brief Derivative of an individual variable in the set.
     */
    template <int i, class Evaluators>
    auto derivative(Evaluators const& evals) const
    {
      int const spaceIndex = Descriptions::template spaceIndex<i>;
      return component<i>(*this).derivative(boost::fusion::at_c<spaceIndex>(evals));
    }
    
    /**
     * \brief Derivative of an individual variable in the set.
     */
    template <int i>
    auto derivative(Cell<GridView> const& cell, LocalPosition<GridView> const& xi) const
    {
      return component<i>(*this).derivative(cell,xi);
    }
    
    /**
     * \brief derivative of an individual variable in the set. 
     * \WARNING This method is comparatively DEAD SLOW.  
     */
    template <int i>
    auto derivative(GlobalPosition<GridView> const& x) const
    {
      return component<i>(*this).derivative(x);
    }

    /**
     * \brief Hessian of an individual variable in the set.
     */
    template <int i, class Evaluators>
    auto hessian(Evaluators const& evals) const
    {
      int const spaceIndex = Descriptions::template spaceIndex<i>;
      return component<i>(*this).hessian(boost::fusion::at_c<spaceIndex>(evals));
    }
    
    /**
     * \brief Hessian of an individual variable in the set.
     */
    template <int i>
    auto hessian(Cell<GridView> const& cell, LocalPosition<GridView> const& xi) const
    {
      return component<i>(*this).hessian(cell,xi);
    }
    
    /**
     * \brief Hessian of an individual variable in the set. 
     * \WARNING This method is comparatively DEAD SLOW.  
     */
    template <int i>
    auto hessian(GlobalPosition<GridView> const& x) const
    {
      return component<i>(*this).hessian(x);
    }
    
    /// @}

    /// Descriptions of variable set, of type VariableSetDescription (lots of useful infos)
    Descriptions const& descriptions;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  /**
   * \ingroup variables
   * \brief A class that stores information about a set of variables.
   *
   * PDE systems usually involve a set of variables \f$ u_1, \dots, u_n \f$, e.g., velocity and 
   * pressure in the Stokes equation. This class describes the properties of such variable sets 
   * and provides meta information such as types of variable sets, their function spaces, or
   * their dimension, as well as convenience methods for constructing variable sets or coefficient
   * vector sets.
   * 
   * 
   *
   *
   * \tparam SpaceList a boost::fusion container of FEFunctionSpace s
   * \tparam VariableList a boost::fusion container of VariableDescription s. Note that
   * the variable ids will be ignored and ids according to the position of the variable in the vector 
   * will be used instead.
   */
  template <class SpaceList, class VariableList>
  class VariableSetDescription
  {
  public:
    /// type of boost::fusion::vector of FEFunctionSpace s needed
    typedef SpaceList Spaces;
    
    /// type of boost::fusion vector of VariableDescription s
    using Variables = typename Variables_Detail::ImplicitIdVariableDescriptions<VariableList>::type;

    /**
     * \brief Type that contains a set of variables inside a boost vector, together with all the data.
     * 
     * This is a boost::fusion::vector with Dune::BlockVector entries.
     */
    typedef typename Variables_Detail::VariableRangeCreator<Variables,ConstructElement<Spaces>>::type RepresentationData;

    /**
     * \brief Type that contains a set of variable values with some
     * functionality, such as simple vector arithmetic and so on
     * (cf. class \ref VariableSet)
     *
     * Construct an object r of this type from a given VariableSetDescription vsd by:
     * \code
     * VariableSetDescription::VariableSet x(vsd);
     * \endcode
     * or by
     * \code
     * auto x = vsd.variableSet();
     * \endcode
     */
    typedef Kaskade::VariableSet<VariableSetDescription<Spaces,VariableList>> VariableSet;
    
    /// scalar field type
    typedef typename SpaceType<Spaces,0>::type::Scalar   Scalar;
    /// Grid type
    typedef typename SpaceType<Spaces,0>::type::Grid     Grid;
    /// Grid view type
    typedef typename SpaceType<Spaces,0>::type::GridView GridView;
    /// IndexSet type
    typedef typename SpaceType<Spaces,0>::type::IndexSet IndexSet;

    /// Number of variables in this set
    static int const noOfVariables = boost::fusion::result_of::size<Variables>::type::value;
    
    /**
     * \brief minimal continuity of any variable in the set 
     * 
     * This is the number of how often all the functions are continuously differentiable. 0 means globally continuous 
     * with discontinuous derivatives. -1 denotes bounded but discontinuous functions.
     */
    static int const continuity = Variables_Detail::Continuity<Variables,Spaces>::value;
    
    /**
     * \brief The boost::fusion sequence of evaluators belonging to this variable set.
     */
    typedef typename boost::fusion::result_of::as_vector<
          typename boost::fusion::result_of::transform<Spaces, GetEvaluators<ShapeFunctionCache<Grid,Scalar>> >::type
          >::type Evaluators;

    /**
     * \name Constructors
     * @{
     */

    /// Constructor, giving each variable a name (e.g. for output purposes)
    template <class RAIter>
    VariableSetDescription(Spaces const& spaces_, RAIter nameIt):
      VariableSetDescription(spaces_)
    {
      std::copy(nameIt,nameIt+noOfVariables,names.begin());
    }

    /**
     * \brief Constructor.
     * \param spaces a heterogeneous sequence of pointers to the (const) involved spaces.
     *               This is copied internally and need not persist.
     * \param names contains a name for each variable.
     */
    VariableSetDescription(Spaces const& spaces_, std::vector<std::string> names_)
    : VariableSetDescription(spaces_)
    {
      assert(names.size()==names_.size());
      std::copy(names_.begin(),names_.end(),names.begin());
    }

    /**
     * \brief Constructor specifying for each variable both a name and a role.
     */
    template <class RAIter>
    VariableSetDescription(Spaces const& spaces_, RAIter nameIt, RAIter roleIt)
    : VariableSetDescription(spaces_,nameIt)
    {
      std::copy(roleIt,roleIt+noOfVariables,roles.begin());
    }

    /**
     * \brief Constructor specifying for each variable both a name and a role.
     */
    VariableSetDescription(Spaces const& spaces_, std::vector<std::string> names_, std::vector<std::string> roles_)
    : VariableSetDescription(spaces_,names_)
    {
      std::copy(roles_.begin(),roles_.end(),roles.begin());
    }

    /// Constructor without naming the variables
    VariableSetDescription(Spaces const& spaces_)
    : spaces(spaces_)
    , gridView(boost::fusion::at_c<0>(spaces_)->gridView())
    , indexSet(boost::fusion::at_c<0>(spaces_)->indexSet())
    {}

    /**
     * @}
     */

    /**
     * \name Size and Dimension queries
     * @{
     */

    /**
     * \brief Computes the total number of scalar degrees of freedom collected in the variables [first,last).
     * 
     * For each affected variable, the number of components is multiplied with the degrees of freedom of its FE space.
     */
    size_t degreesOfFreedom(int first=0, int last=noOfVariables) const
    {
      return degreesOfFreedom(spaces,first,last);
    }

    /**
     * \brief Computes the total number of scalar degrees of freedom collected in the variables [first,last).
     * 
     * For each affected variable, the number of components is multiplied with the degrees of freedom of its FE space.
     */
    static size_t degreesOfFreedom(Spaces const& spaces, int first, int last)
    {
      auto dimensions = variableDimensions(spaces);
      assert(0<=first && first<=last && last<=dimensions.size());
      return std::accumulate(begin(dimensions)+first,begin(dimensions)+last,0);
    }      
    
    /**
     * \brief Computes for each variable the number of scalar degrees of freedom.
     */
    static std::array<size_t,noOfVariables> variableDimensions(Spaces const& spaces, int first=0, int last=noOfVariables)
    {
      std::array<size_t,noOfVariables> dimensions;
      boost::fusion::for_each(Variables(),ComputeDimension<noOfVariables>(dimensions,spaces));
      return dimensions;
    }

    /**
     * \brief number of components of i'th variable
     */
    template <int i>
    constexpr static int components = boost::fusion::result_of::value_at_c<Variables,i>::type::m;

    /// DEPRECATED, use components instead
    template <int idx>
    struct [[deprecated]] Components { static int const m = boost::fusion::result_of::value_at_c<Variables,idx>::type::m; };

    /**
     * @}
     */

    /**
     * \brief space index of the i'th variable
     *
     * For querying the space index of variables, a more convenient method is to use the global spaceIndex
     * templated value:
     * \code
     * int varSpaceIdx = spaceIndex<VariableSetDesc,varIdx>;
     * \endcode
     */
    template <int idx>
    constexpr static int spaceIndex = boost::fusion::result_of::value_at_c<Variables,idx>::type::spaceIndex;

    /**
     * \brief DEPRECATED
     */
    template <int idx>
    struct [[deprecated]] SpaceIndex { static int const spaceIndex = boost::fusion::result_of::value_at_c<Variables,idx>::type::spaceIndex; };

    
    /**
     * \brief The type of finite element space of the i-th variable.
     */
    template <int idx>
    using Space = std::remove_pointer_t<typename boost::fusion::result_of::value_at_c<Spaces,SpaceIndex<idx>::spaceIndex>::type>;
    
    /**
     * \name Creation of variable sets and coefficient vectors
     * @{
     */

    /**
     * \brief Typedefs for coefficient vector representations of contiguous subranges of the variables.
     * 
     * Note that the representation defined by this class does not consist of FE functions.
     *
     * Since in contrast to FE functions the linear algebra coefficient
     * vectors are invalidated during mesh adaptation, the coefficient
     * vectors should be copied to FE functions before mesh refinement.
     */
    template <int first=0, int last=noOfVariables>
    class CoefficientVectorRepresentation 
    {
      typedef Variables_Detail::VariableRangeCreator<Variables,ConstructCoefficientVector<Spaces>,first,last> Creator;

    public:
      /**
       * \brief A linear product space of FE coefficient vectors for the variables.
       * Note that this is not a linear space of FE functions, only of coefficient vectors. 
       * In contrast to FE functions, coefficient vectors are invalidated during mesh
       * adaptation.
       */
      typedef LinearProductSpace<Scalar,typename Creator::type> type;
      
      /**
       * \brief Returns a (small) initializer object for coefficient vectors.
       * 
       * The coeffcient vectors initialized with the returned initializer object are set to zero.
       */
      static typename Creator::View init(Spaces const& spaces) 
      {
        return Creator::apply(ConstructCoefficientVector<Spaces>(spaces));
      }
    };
    
    /**
     * \brief Returns a zero initialized coefficient vector for variables [first,last[.
     */
    template <int first=0, int last=noOfVariables>
    auto zeroCoefficientVector() const 
    {
      using CVec = CoefficientVectorRepresentation<first,last>;
      return typename CVec::type(CVec::init(spaces));
    }
    
    /**
     * \brief Returns a zero-initialized variable set.
     */
    VariableSet variableSet() const 
    {
      return VariableSet(*this);
    }

    /**
     * @}
     */

    /**
     * \brief A heterogeneous sequence of pointers to (const) spaces involved.
     */
    Spaces             spaces;
    
    /**
     * \brief The grid view on which the variables live.
     */
    GridView const&    gridView;
    
    /**
     * \brief The gridManager keeping the grid on which we operate.
     */
    GridManagerBase<Grid> const& gridManager() const { return boost::fusion::at_c<0>(spaces)->gridManager(); }

    /// index set
    IndexSet const&    indexSet;
    
    /**
     * \brief A sequence of variable names
     */
    std::array<std::string,noOfVariables>  names;

    std::array<std::string,noOfVariables>  roles;

  private:

    struct CheckHasEqualGridAndIndexSet
    {
      CheckHasEqualGridAndIndexSet(GridManagerBase<Grid> const* gp_, IndexSet const* isp_):
        gp(gp_), isp(isp_)
      {}

      template <class Space>
      void operator()(Space const* space) const {
        assert(&(space->gridManager())== gp);
        assert(&(space->indexSet()) == isp);
      }

      GridManagerBase<Grid> const* gp;
      IndexSet const* isp;
    };

    template <int nvars>
    struct ComputeDimension
    {
      ComputeDimension(std::array<size_t,nvars>& dims, Spaces const& spaces_):
        dims_(dims), spaces(spaces_)
      {}

      template <class Variable>
      void operator()(Variable const& /* v */) const
      {
        assert(dims_.size()>Variable::id);
        dims_[Variable::id] = Variable::m * boost::fusion::at_c<Variable::spaceIndex>(spaces)->degreesOfFreedom();
      }

      std::array<size_t,nvars>& dims_;
      Spaces const&             spaces;
    };

  };

  /**
   * \ingroup variables
   * \brief The index of the finite element space of the variable with given index within the list of spaces.
   * \relates VariableSetDescription
   * \relates VariableSetDescription::SpaceIndex
   */
  template <class VarSetDesc, int varIdx>
  constexpr int spaceIndex = VarSetDesc::template spaceIndex<varIdx>;
  
  //---------------------------------------------------------------------
  
  
  /**
   * \ingroup variables 
   * \brief Creates a list of finite element spaces in the required format.
   * Usage example (for Stokes flow):
   * \code
   * auto spaces = makeSpaceList(&velocitySpace,&pressureSpace);
   * \endcode
   * This creates a boost::fusion vector of pointers to const spaces, as expected by
   * \ref VariableSetDescription.
   *
   * \relates VariableSetDescription
   */
  template <class ...Spaces>
  auto makeSpaceList(Spaces*... spaces)
  {
    return boost::fusion::make_vector(const_cast<Spaces const*>(spaces)...);
  }
  
  /**
   * \ingroup variables 
   * \brief Creates a variable set description from given spaces and variables.
   * \param spaces a boost::fusion vector of pointers to the (const) spaces referenced by the variables 
   * \param vars a boost::fusion sequence of \ref Variable objects describing the variables in the set
   * 
   * \code
   * auto variableSetDescription = makeVariableSetDescription(makeSpaceList(&velocitySpace,&pressureSpace),
   *                                                          boost::fusion::make_vector(Variable<SpaceIndex<0>,Components<2>>("u"),
   *                                                                                     Variable<SpaceIndex<1>,Components<1>>("p")));
   * \endcode
   * \relates VariableSetDescription
   */
  template <class SpaceList, class VariableList>
  VariableSetDescription<SpaceList,VariableList> makeVariableSetDescription(SpaceList const& spaces,
                                                                            VariableList const& vars)
  {
    using namespace boost::fusion;
    
    // Extract the names of the variables.
    std::vector<std::string> names;
    for_each(vars,[&](auto var) { names.push_back(var.name); });
    // If (some) variable names are not specified, assign default ones.
    for (int i=0; i<names.size(); ++i)
      if (names[i].empty())
        names[i] = "var" + paddedString(i,2);
      
    return VariableSetDescription<SpaceList,VariableList>(spaces,names);
  }
  

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  

  /**
   * \ingroup variables
   * \brief A function evaulating all functions in a variable set using the provided method.
   * 
   * \tparam VariableSet a heterogeneous sequence of variable descriptions (defines the mapping of functions to evaluators)
   * \tparam Functions a heterogeneous sequence of FE functions (or function views)
   * \tparam Evaluators the heterogeneous container of evaluators for the variable set 
   * \tparam Method a polymorphic Functor type with call operator (Function const&, Evaluator const&)
   * 
   * \param vars the variable set to be evaluated.
   * \param eval the evaluators to be used for function evaluation
   * \param method the functor performing the actual evaluation
   * 
   *  evaluateVariables valueMethod ValueMethod derivativeMethod DerivativeMethod
   */
  template <class VariableDescriptions, class Functions, class Evaluators, class Method>
  auto evaluateFunctions(Functions const& fs,
                         Evaluators const& eval,
                         Method const& method = Method())
  {
    using namespace boost::fusion;
    
    // Define the functor that evaluates the individual variables.
    auto doEval = [&] (auto const& pair)
    {
      // Extract the space index from the second component of the pair provided by zip.
      static int const sid = std::remove_reference_t<typename result_of::value_at_c<std::remove_reference_t<decltype(pair)>,1>::type>::spaceIndex;
      // Evaluate with evaluator corresponding to the function's space.
      return method(at_c<0>(pair),at_c<sid>(eval));
    };
    
    // We need this as vector type to default construct an object for transform. This way we can submit
    // views (such as joint_view) as template parameter.
    using VarDesc = typename result_of::as_vector<VariableDescriptions>::type;
    
    // For evaluation, the corresponding evaluator must be used based on the space index of the variables. To provide this,
    // we zip the functions with their variable descriptions to pairs from which we can extract both function and space index.
    return as_vector(transform(zip(fs,VarDesc()),doEval));
  }
  
  /**
   * \ingroup variables
   * \brief A function evaulating all functions in a variable set using the provided method.
   * 
   * \tparam VariableSet the type of the variable set to be evaluated (i.e. VariableSet<...>)
   * \tparam Method a polymorphic Functor type with call operator (Function const&, Evaluator const&)
   * 
   * \param vars the variable set to be evaluated.
   * \param eval a heterogeneous sequence of evaluators to be used for function evaluation
   * \param method the functor performing the actual evaluation (usually valueMethod or derivativeMethod)
   *
   * \return a boost::fusion::vector containing the results of variable evaluation
   * 
   * \see EvaluateVariables, evaluateFunctions, valueMethod, ValueMethod, derivativeMethod, DerivativeMethod
   */
  template <class VariableSet, class Method>
  auto evaluateVariables(VariableSet const& vars,
                         typename VariableSet::Descriptions::Evaluators const& eval,
                         Method const& method = Method())
  {
    return evaluateFunctions<typename VariableSet::Descriptions::Variables>(vars.data,eval,method);
  }
  
  /**
   * \ingroup variables
   * \brief The type of evaluated variables as returned by \ref evaluateVariables.
   * 
   * This is a boost::fusion sequence of functions' values resulting from evaluating all functions in a given variable set.
   * \tparam VariableSet the type of the variable set (i.e. VariableSet<...>)
   * \tparam Method a polymorphic Functor type with call operator (Function const&, Evaluator const&) used for evaluating.
   * 
   * \related evaluateVariables
   */
  template <class VariableSet, class Method>
  using EvaluateVariables = decltype(evaluateVariables(std::declval<VariableSet>(),
                                                       std::declval<typename VariableSet::Descriptions::Evaluators>(),
                                                       std::declval<Method>()));


  /**
   * \ingroup variables
   * \brief Helper method for evaluating whole variable sets.
   * 
   * Use this in calling evaluateVariables when evaluating the functions' values.
   * \relates evaluateVariables
   */
  constexpr auto valueMethod = [](auto const& f, auto const& eval) { return f.value(eval); };
  
  /**
   * \ingroup variables
   * \brief Helper method for evaluating whole variable sets.
   * 
   * Use this in calling evaluateVariables when evaluating the functions' derivatives.
   * \relates evaluateVariables
   */
  constexpr auto derivativeMethod = [](auto const& f, auto const& eval) { return f.derivative(eval); };
  
  /**
   * \ingroup variables
   * \brief The type of the valueMethod helper functor.
   */
  using ValueMethod = decltype(valueMethod);
  
  /**
   * \ingroup variables
   * \brief The type of the derivativeMethod helper functor.
   */
  using DerivativeMethod = decltype(derivativeMethod);

} // end of namespace Kaskade




//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------

#endif
