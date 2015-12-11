/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef VIEWER_H
#define VIEWER_H

#include <QGLWidget>

class QMouseEvent;

/*!
 * \brief The Viewer class
 *
 *
 *
 * \see Trackball implementation from http://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php
 */
class Viewer : public QGLWidget
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();

    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);

private:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);

    void draw();

    int _xRot, _yRot, _zRot;

    QPoint _lastPos;
};

#endif // VIEWER_H
