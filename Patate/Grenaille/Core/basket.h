#ifndef _GRENAILLE_BASKET_
#define _GRENAILLE_BASKET_


namespace Grenaille{




  template <class, class, typename T> class Forward: public T {};


#define BASKET_TP(I) template <class, class, typename> class Ext##I = Forward
  /*!

    \brief 
    \todo Comment

   */
  template < class P, class W, template <class, class, typename>class Fit, 
	     BASKET_TP(0),  BASKET_TP(1), BASKET_TP(2), BASKET_TP(3), BASKET_TP(4), BASKET_TP(5) >
  class Basket
    : public Ext5<P,W,Ext4<P,W, Ext3<P,W, Ext2<P,W, Ext1<P,W, Ext0<P,W, Fit<P,W,void> > > > > > >
  {
    //  private:
    //    typedef Ext5<S,D,W,Ext4<S,D,W, Ext3<S,D,W, Ext2<S,D,W, Ext1<S,D,W, Ext0<S,D,W, Fit<S,D,W> > > > > > > Base;
  }; // class Basket

#undef BASKET_TP

}// namespace Grenaille


#endif
