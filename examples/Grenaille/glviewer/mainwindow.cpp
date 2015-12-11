/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/


#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>

#include "Patate/common/surface_mesh/objReader.h"

typedef Viewer::Mesh Mesh;


bool noopErrorCallback(const std::string& /*msg*/, unsigned /*line*/, void* /*ptr*/) {
    return false;
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionOpen_File_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(
                this,
                "Load obj file",
                ".",
                QString("3d Mesh (*.obj)"));

    if (! fileName.isEmpty()) {
        Mesh mesh;
        {
            std::ifstream in(fileName.toStdString().c_str());
            PatateCommon::OBJReader<Mesh> reader;
            reader.setErrorCallback(PatateCommon::defaultErrorCallback,
                                    noopErrorCallback, NULL);
            reader.read(in, mesh);
        }
        ui->_viewer->setMesh(mesh);
    }
}
