/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef VIEWER_H
#define VIEWER_H

#ifdef __APPLE__
    #include <OpenGL/gl3.h>
    #include <OpenGL/gl3ext.h>

    #define __gl_h_

#elif defined _MSC_VER
    #include <Windows.h>
    #include <GL/gl.h>
#else
    #define GL_GLEXT_PROTOTYPES

    #include <GL/gl.h>
    #include <GL/glext.h>
#endif

#include <QGLWidget>
//#include <QGLFunctions>


#include "Patate/common/gl_utils/glmesh.h"
#include <QOpenGLShaderProgram>

class QMouseEvent;

/*!
 * \brief The Viewer class
 *
 *
 *
 * \see Trackball implementation from http://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php
 */
class Viewer : public QGLWidget/*, protected QGLFunctions*/
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();

    typedef PatateCommon::GLTri3DMesh Mesh;

    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseDoubleClickEvent(QMouseEvent *event);

    inline void setMesh(Mesh* mesh) { _mesh = mesh; update(); }

public slots:
    void setScale(double scale) { _scale = scale; update(); }

private:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);

    //! \brief Enable/Disable auto refresh
    void triggerAutoRefresh(bool status);

    void draw(Mesh *mesh);
    void drawPicked();
    void prepareShaders();
    void prepareFBO(int w, int h);
    void updateTransformationMatrix(bool reset = true);

    int _xRot, _yRot, _zRot, _zoom;

    QPoint _lastPos;
    QOpenGLShaderProgram _program, _pickProgram;
    struct {
        int vertex;
        int normal;
        int transform;
        int lightPos;
    } _progLocation;
    struct {
        int vertex;
        int normal;
        int ids;
        int transform;
    } _pickingProgLocation;
    GLuint _pickingFBOLocation, _pickingDepth, _pickingTexture;
    QVector3D _lightPos;
    QMatrix4x4 _transform;
    bool _programInitialized;

    int _pickedPointId;
    typename Mesh::Vector _pickedPoint;

    double _scale;

    Mesh *_mesh;
    Mesh _unitSphere;
    QTimer *_refreshTimer; // automatically update the view
};

#endif // VIEWER_H
