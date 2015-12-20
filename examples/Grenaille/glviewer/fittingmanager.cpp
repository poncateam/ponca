#include "fittingmanager.h"

FittingManager::FittingManager(QObject *parent) :
    QObject(parent),
    _fitType(FittingManager::UNSUPPORTED),
    _mesh(NULL)
{
}
void
FittingManager::setBasketType(FIT_TYPE type){
    _fitType = type;
}
