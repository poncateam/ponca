/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/

#include <QApplication>
#include <QDesktopWidget>
#include "mainwindow.h"

int main(int argc, char ** argv)
{
    QApplication app(argc, argv);
    MainWindow window;
    window.show();

    return app.exec();
}
