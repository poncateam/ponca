/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#include "viewer.h"
#include <QMouseEvent>
#include <QOpenGLShaderProgram>

/*!
 * \brief Viewer::Viewer
 * \param parent
 */
Viewer::Viewer(QWidget *parent) :
    QGLWidget(parent),
    _xRot(0),
    _yRot(0),
    _zRot(0),
    _zoom(10),
    _lastPos(QPoint(0,0)),
    _programInitialized(false)
{
}

Viewer::~Viewer()
{
}


static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360)
        angle -= 360 * 16;
}

void Viewer::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != _xRot) {
        _xRot = angle;
        update();
    }
}

void Viewer::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != _yRot) {
        _yRot = angle;
        update();
    }
}

void Viewer::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != _zRot) {
        _zRot = angle;
        update();
    }
}


void Viewer::initializeGL()
{
    qglClearColor(Qt::white);

    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);

//    _lightPos.setX(0);
//    _lightPos.setY(0);
//    _lightPos.setZ(1);

    prepareShaders();
}

void Viewer::prepareShaders() {
    _program.addShaderFromSourceFile(QOpenGLShader::Vertex, "default.vert");
    _program.addShaderFromSourceFile(QOpenGLShader::Fragment, "default.frag");

    _program.link();
    _program.bind();

    _progLocation.vertex    = _program.attributeLocation("vertex");
    _progLocation.normal    = _program.attributeLocation("normal");
    _progLocation.transform = _program.uniformLocation("transform");
    //_progLocation.lightPos  = _program.uniformLocation("lightPos");

    _program.release();
    _programInitialized = true;
}


void Viewer::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    _transform.setToIdentity();
    _transform.scale(_zoom * 0.1);
    //_transform.translate(0.0, 0.0, 0.0);
    _transform.rotate( -_xRot / 16.0, 1.0, 0.0, 0.0);
    _transform.rotate( -_yRot / 16.0, 0.0, 1.0, 0.0);
    _transform.rotate( _zRot / 16.0, 0.0, 0.0, 1.0);
    draw();
}

void Viewer::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(-2, +2, -2, +2, 0.1, 100.0);
#else
    glOrtho(-2, +2, -2, +2, 0.01, 200.0);
#endif
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    _lastPos = event->pos();
}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - _lastPos.x();
    int dy = event->y() - _lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(_xRot + 8 * dy);
        setYRotation(_yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(_xRot + 8 * dy);
        setZRotation(_zRot + 8 * dx);
    }

    _lastPos = event->pos();
}

void Viewer::wheelEvent(QWheelEvent * event){
    static const int unit = 120;
    _zoom = std::max(1, _zoom + event->angleDelta().y()/unit );
    update();
}

void Viewer::draw(){
     // draw object
    _program.bind();
    _program.enableAttributeArray(_progLocation.vertex);
    _program.enableAttributeArray(_progLocation.normal);
    _program.setUniformValue("transform", _transform);
    //_program.setUniformValue("lightPos", _lightPos);
    _mesh.draw();
    _program.disableAttributeArray(_progLocation.normal);
    _program.disableAttributeArray(_progLocation.vertex);
    _program.release();
}
