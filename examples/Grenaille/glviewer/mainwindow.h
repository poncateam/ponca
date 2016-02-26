/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "fittingmanager.h"

#include <QMainWindow>

namespace Ui {
class MainWindow;
}
namespace PatateCommon {
class GLTri3DMesh;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_actionOpen_File_triggered();
    void on__paramBasketType_currentIndexChanged(int index);

    void updateNeighborhoodMesh();

private:
    Ui::MainWindow *ui;
    FittingManager _manager;
    PatateCommon::GLTri3DMesh* _mesh;
};

#endif // MAINWINDOW_H
