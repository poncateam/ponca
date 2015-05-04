/*
Copyright (C) 2013 Nicolas Mellado <nmellado0@gmail.com>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef CIRCLEGL_H
#define CIRCLEGL_H

#include "GL/glew.h"
#include "Eigen/Core"

namespace internal
{

/*!
* Internal class providing static assertion.
* Designed using the static assertion mechanism of Eigen (http://eigen.tuxfamily.org/)
* \FIXME Visual Studio compatibility must be taken into account. Is maybe an occasion to put static assert in patate ?
*/
template<bool condition>
struct StaticCircleAssertion{};

template<>
struct StaticCircleAssertion<false>
{
	enum
	{
		CIRCLEGL_IS_DEFINED_ONLY_FOR_BASKET_OPERATING_ON_2D_ELEMENTS
	};
};

}

/*!
* Display an algebraic Sphere in dimension 2 as a circle using OpenGL
*/
template <class DataPoint, class _WFunctor, typename T>
class CircleGL : public T
{
private:
    typedef T Base;
    GLuint m_vbo [1];
    unsigned int m_slices;

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE,
        PROVIDES_GL_DISPLAY
    };

public:
    inline void initGL(bool _revertY = false, unsigned int _slices = 32);
    inline void displayGL() const;

}; //class CircleGL





////////////////////////////////////////////////
/// Implementation
///

template <class DataPoint, class _WFunctor, typename T>
void
CircleGL<DataPoint, _WFunctor, T>::initGL(bool _revertY, unsigned int _slices)
{
    if (internal::StaticCircleAssertion<(DataPoint::Dim != 2)>::CIRCLEGL_IS_DEFINED_ONLY_FOR_BASKET_OPERATING_ON_2D_ELEMENTS)
    {}

    // setup and cleaning
    m_slices = _slices;
    glDeleteBuffers(1, m_vbo);


    // Create Circle VBO arrays on CPU
    Eigen::Vector2f data [m_slices];
    const double step   =  2.f*M_PI / (m_slices-1);
    const Eigen::Vector2f center = Base::center();
    const double radius           = Base::radius();

    for (int i = 0; i != m_slices; i++)
    {
        const double angle = i * step;

        Eigen::Vector2f &current = data[i];

        current << radius * cos(angle), radius * sin(angle);
        current += center;

        if (_revertY)
        {
            current(1) = 1.f - current(1);
        }
    }


    // Transfert VBO to GPU
    glGenBuffers(1, m_vbo);

    // load coordinate
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
    glBufferData(GL_ARRAY_BUFFER,
                m_slices*sizeof(Eigen::Vector2f),
                data[0].data(),
                GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER,0);
}

template <class DataPoint, class _WFunctor, typename T>
void
    CircleGL<DataPoint, _WFunctor, T>::displayGL() const
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[0]);
    glVertexPointer( 2, GL_FLOAT, 0, 0);
    glDrawArrays(GL_LINE_STRIP, 0, m_slices);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glDisableClientState(GL_VERTEX_ARRAY);
}


#endif // CIRCLEGL_H
