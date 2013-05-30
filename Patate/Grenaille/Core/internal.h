/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_INTERNAL_
#define _GRENAILLE_INTERNAL_


namespace Grenaille
{
  namespace internal{

    /*!
     \class FixedSizeArray
     \brief Internal class to store data in an array with known size at compile time.
    */
    template <typename _DataType, int _Size>
    class FixedSizeArray{
    private:
      _DataType _data [_Size];
      
    public:
      MULTIARCH static inline int size() {return _Size;}
      
      MULTIARCH inline       _DataType& operator[] (int id)       {return _data[id];}
      MULTIARCH inline const _DataType& operator[] (int id) const {return _data[id];}
      

      // use template function
      template <class Scalar> MULTIARCH
      FixedSizeArray<_DataType,_Size> operator*(const Scalar& s){
        FixedSizeArray<_DataType,_Size> result;
        for(unsigned int i = 0; i != _Size; i++) result[i] = s*_data[i];
        return result;
      }
      
      MULTIARCH inline _DataType* data() {return _data;}         
    };  //class FixedSizeArray
  }  //namespace internal
}  //namespace Grenaille


#endif
