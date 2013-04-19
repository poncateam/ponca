#ifndef _GRENAILLE_BASKET_
#define _GRENAILLE_BASKET_


namespace Grenaille{



  namespace internal{
    template <class, class, typename T> class Forward: public T {};
  }


#define BASKET_TP(I) template <class, class, typename> class Ext##I = internal::Forward
  /*!

    \brief 
    \todo Comment

   */
  template < class P, class W, template <class, class, typename>class Fit, 
	     BASKET_TP(0), BASKET_TP(1), BASKET_TP(2), BASKET_TP(3), BASKET_TP(4), BASKET_TP(5), BASKET_TP(6), BASKET_TP(7), BASKET_TP(8), BASKET_TP(9), BASKET_TP(10), BASKET_TP(11) >
  class Basket
    : public Ext11<P,W, Ext10<P,W, Ext9<P,W, Ext8<P,W, Ext7<P,W, Ext6<P,W, Ext5<P,W, Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, Fit<P,W,void> > > > > > > > > > > > > 
  {
  }; // class Basket

#undef BASKET_TP

}// namespace Grenaille


#endif
