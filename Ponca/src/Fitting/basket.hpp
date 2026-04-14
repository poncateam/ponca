namespace Ponca
{
    /*
        Keep the following functions in a separate file. In order for functions to 
        be precompiled by Ponca, we need them not to be inline. 
        In particular, their body can not be inside the class definitions.
    */

// Shortcut for functions needing this type
#define DiffArray typename internal::BasketDiffAggregate<BasketType, Type, Ext0, Exts...>::type::ScalarArray

    template <typename BasketType, int Type, template <class, class, int, typename> class Ext0,
              template <class, class, int, typename> class... Exts>
    PONCA_MULTIARCH void BasketDiff<BasketType, Type, Ext0, Exts...>::init()
    {
        Base::init();
    }

    template <typename BasketType, int Type, template <class, class, int, typename> class Ext0,
              template <class, class, int, typename> class... Exts>
    PONCA_MULTIARCH void BasketDiff<BasketType, Type, Ext0, Exts...>::startNewPass()
    {
        Base::startNewPass();
    }

    template <typename BasketType, int Type, template <class, class, int, typename> class Ext0,
              template <class, class, int, typename> class... Exts>
    PONCA_MULTIARCH void BasketDiff<BasketType, Type, Ext0, Exts...>::finalize()
    {
        Base::finalize();
    }

    template <typename BasketType, int Type, template <class, class, int, typename> class Ext0,
              template <class, class, int, typename> class... Exts>
    PONCA_MULTIARCH bool BasketDiff<BasketType, Type, Ext0, Exts...>::addNeighbor(const BSKP& _nei)
    {
        // compute weight
        auto neiFilterOutput = Base::getNeighborFilter()(_nei);
        typename Base::ScalarArray dw;

        if (neiFilterOutput.first > Scalar(0.))
        {
            Base::addLocalNeighbor(neiFilterOutput.first, neiFilterOutput.second, _nei, dw);
            return true;
        }
        return false;
    }

    template <typename BasketType, int Type, template <class, class, int, typename> class Ext0,
              template <class, class, int, typename> class... Exts>
    PONCA_MULTIARCH void BasketDiff<BasketType, Type, Ext0, Exts...>::addLocalNeighbor(BSKS _w, const BSKV& _localQ, const BSKP& _nei, DiffArray& _dw)
    {
        Base::addLocalNeighbor(_w, _localQ, _nei, _dw);
    }

    template <class P, class NF, template <class, class, typename> class Ext0,
              template <class, class, typename> class... Exts>
    PONCA_MULTIARCH void Basket<P, NF, Ext0, Exts...>::init()
    {
        Base::init();
    }

    template <class P, class NF, template <class, class, typename> class Ext0,
              template <class, class, typename> class... Exts>
    PONCA_MULTIARCH void Basket<P, NF, Ext0, Exts...>::startNewPass()
    {
        Base::startNewPass();
    }

    template <class P, class NF, template <class, class, typename> class Ext0,
              template <class, class, typename> class... Exts>
    PONCA_MULTIARCH void Basket<P, NF, Ext0, Exts...>::finalize()
    {
        Base::finalize();
    }

    template <class P, class NF, template <class, class, typename> class Ext0,
              template <class, class, typename> class... Exts>
    PONCA_MULTIARCH bool Basket<P, NF, Ext0, Exts...>::addNeighbor(const P& _nei)
    {
        // compute weight
        auto neiFilterOutput = Base::getNeighborFilter()(_nei);

        if (neiFilterOutput.first > Scalar(0.))
        {
            Base::addLocalNeighbor(neiFilterOutput.first, neiFilterOutput.second, _nei);
            return true;
        }
        return false;
    }

    template <class P, class NF, template <class, class, typename> class Ext0,
              template <class, class, typename> class... Exts>
    PONCA_MULTIARCH void Basket<P, NF, Ext0, Exts...>::addLocalNeighbor(typename P::Scalar _w, const typename P::VectorType& _localQ, const P& _nei)
    {
        Base::addLocalNeighbor(_w, _localQ, _nei);
    }

#undef DiffArray
}