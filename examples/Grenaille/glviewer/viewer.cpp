/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#include "viewer.h"
#include "Patate/common/surface_mesh/objReader.h"


#include <istream>

#include <QTimer>
#include <QMouseEvent>
#include <QOpenGLShaderProgram>

#define ASSET_SPHERE_MESH_FILENAME "unitSphere.obj"


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
    _programInitialized(false),
    _pickedPointId(-1),
    _scale(0.02),
    _mesh(NULL),
    _refreshTimer(new QTimer(this))
{
    _refreshTimer->setInterval(25);
    _refreshTimer->setSingleShot(false);
    connect(_refreshTimer, SIGNAL(timeout()), this, SLOT(update()));
    triggerAutoRefresh(false);

    // load unit sphere mesh
    std::ifstream in(ASSET_SPHERE_MESH_FILENAME);
    PatateCommon::OBJReader<Mesh> reader;
    reader.read(in, _unitSphere);

    updateTransformationMatrix();
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
        updateTransformationMatrix();
    }
}

void Viewer::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != _yRot) {
        _yRot = angle;
        updateTransformationMatrix();
    }
}

void Viewer::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != _zRot) {
        _zRot = angle;
        updateTransformationMatrix();
    }
}

void Viewer::triggerAutoRefresh(bool status) {
    if (status){
        _refreshTimer->start();
        update();
    }
    else
        _refreshTimer->stop();
}


void Viewer::initializeGL()
{

    glEnable(GL_DEPTH_TEST);

    prepareShaders();

    glGenFramebuffers(1, &_pickingFBOLocation);
    glGenTextures(1, &_pickingTexture);

    updateProjectionMatrix();
}

void Viewer::prepareShaders() {
    _program.addShaderFromSourceFile(QOpenGLShader::Vertex, "default.vert");
    _program.addShaderFromSourceFile(QOpenGLShader::Fragment, "default.frag");

    _program.link();
    _program.bind();

    _progLocation.vertex      = _program.attributeLocation("vertex");
    _progLocation.normal      = _program.attributeLocation("normal");
    _progLocation.transform   = _program.uniformLocation("transform");
    _progLocation.projection  = _program.uniformLocation("projection");
    //_progLocation.lightPos  = _program.uniformLocation("lightPos");

    _program.release();

    _pickProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "picking.vert");
    _pickProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "picking.frag");

    _pickProgram.link();
    _pickProgram.bind();

    _pickingProgLocation.vertex     = _pickProgram.attributeLocation("vertex");
    _pickingProgLocation.normal     = _pickProgram.attributeLocation("normal");
    _pickingProgLocation.ids        = _pickProgram.attributeLocation("ids");
    _pickingProgLocation.transform  = _pickProgram.uniformLocation("transform");
    _pickingProgLocation.projection = _pickProgram.uniformLocation("projection");

    _pickProgram.release();

    _programInitialized = true;
}

void Viewer::prepareFBO(int w, int h) {
    glBindFramebuffer(GL_FRAMEBUFFER, _pickingFBOLocation);

    // render output
    glBindTexture(GL_TEXTURE_2D, _pickingTexture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32UI, w, h, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, _pickingTexture, 0);

    // depth layer
    glGenRenderbuffers(1, &_pickingDepth);
    glBindRenderbuffer(GL_RENDERBUFFER, _pickingDepth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, _pickingDepth);


    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        std::cerr << "Incomplete FBO Attachement" << std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

}


void Viewer::updateTransformationMatrix(bool reset){
    if(reset) _transform.setToIdentity();

    _transform.scale(_zoom * 0.1);
    _transform.rotate( -_xRot / 16.0, 1.0, 0.0, 0.0);
    _transform.rotate( -_yRot / 16.0, 0.0, 1.0, 0.0);
    _transform.rotate( _zRot / 16.0, 0.0, 0.0, 1.0);
}

void Viewer::updateProjectionMatrix() {
    _projection.setToIdentity();

    const QVector3D s (1., 1., 1./std::sqrt(3.*_zoom));
    _projection.scale(s);

}

void Viewer::paintGL()
{
    qglClearColor(Qt::white);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (_mesh != NULL){
        draw(_mesh);
        if(_pickedPointId != -1)  drawPicked();
    }
}

void Viewer::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    _lastPos = event->pos();
    triggerAutoRefresh(true);
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
    updateTransformationMatrix(true);
    updateProjectionMatrix();
    update();
}

void Viewer::mouseReleaseEvent(QMouseEvent */*event*/)
{
    triggerAutoRefresh(false);
}

void Viewer::mouseDoubleClickEvent(QMouseEvent *event)
{
    _pickedPointId = -1;

    if(_mesh != NULL){
        int w = width();
        int h = height();
        prepareFBO(w,h);

        // bind shader
        _pickProgram.bind();
        _pickProgram.enableAttributeArray(_pickingProgLocation.vertex);
        _pickProgram.enableAttributeArray(_pickingProgLocation.normal);
        _pickProgram.enableAttributeArray(_pickingProgLocation.ids);
        _pickProgram.setUniformValue("transform",  _transform);
        _pickProgram.setUniformValue("projection", _projection);

        // bind fbo
        glBindFramebuffer(GL_FRAMEBUFFER, _pickingFBOLocation);
        glBindTexture(GL_TEXTURE_2D, _pickingTexture);
        glBindRenderbuffer(GL_RENDERBUFFER, _pickingDepth);

        GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
        glDrawBuffers(_pickingFBOLocation, DrawBuffers);
        qglClearColor(Qt::black);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // draw
        _mesh->drawIds();

        // unbind fbo
        glFlush(); // Flush after just in case
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindTexture(GL_TEXTURE_2D, 0);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);

        // unbind shader
        _pickProgram.disableAttributeArray(_pickingProgLocation.vertex);
        _pickProgram.disableAttributeArray(_pickingProgLocation.normal);
        _pickProgram.disableAttributeArray(_pickingProgLocation.ids);
        _pickProgram.release();

        // get texture back
        // SHOULD BE
        // unsigned int* ids = new unsigned int[ w*h ];
        // But it doesn't work as it should...
        // \helpwanted
        float* ids = new float[ w*h ];
        glBindTexture(GL_TEXTURE_2D, _pickingTexture);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, ids);
        glBindTexture(GL_TEXTURE_2D, 0);

        QPoint p = event->pos();
        unsigned int id = ids[p.x() + (h-p.y())*w];
        delete [](ids);
        if (id != 0) {
            _pickedPointId = id;
            _pickedPoint   = _mesh->getVertexMap(id);
            emit selectedPointChanged(_pickedPoint);
        }
        update();
    }
}


void Viewer::draw(Mesh* mesh){
     // draw object
    _program.bind();
    _program.enableAttributeArray(_progLocation.vertex);
    _program.enableAttributeArray(_progLocation.normal);
    _program.setUniformValue("transform", _transform);
    _program.setUniformValue("projection", _projection);
    //_program.setUniformValue("lightPos", _lightPos);
    mesh->draw();
    _program.disableAttributeArray(_progLocation.normal);
    _program.disableAttributeArray(_progLocation.vertex);
    _program.release();
}

void Viewer::drawPicked(){

    _transform.translate(_pickedPoint(0), _pickedPoint(1), _pickedPoint(2));
    _transform.scale(_scale);

    draw(&_unitSphere);

    updateTransformationMatrix();

}
