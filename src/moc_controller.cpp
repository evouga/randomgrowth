/****************************************************************************
** Meta object code from reading C++ file 'controller.h'
**
** Created: Sun Apr 5 17:59:42 2015
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "controller.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'controller.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Controller[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      21,   12,   11,   11, 0x0a,
      44,   12,   11,   11, 0x0a,
      74,   67,   11,   11, 0x0a,
     110,   11,   11,   11, 0x0a,
     117,   11,   11,   11, 0x0a,
     132,   11,   11,   11, 0x0a,
     143,   11,   11,   11, 0x0a,
     154,   11,   11,   11, 0x0a,
     169,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Controller[] = {
    "Controller\0\0filename\0exportOBJ(std::string)\0"
    "importOBJ(std::string)\0params\0"
    "updateParameters(ProblemParameters)\0"
    "quit()\0centerCamera()\0updateGL()\0"
    "makeCone()\0makeFlatCone()\0crush()\0"
};

void Controller::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Controller *_t = static_cast<Controller *>(_o);
        switch (_id) {
        case 0: _t->exportOBJ((*reinterpret_cast< std::string(*)>(_a[1]))); break;
        case 1: _t->importOBJ((*reinterpret_cast< std::string(*)>(_a[1]))); break;
        case 2: _t->updateParameters((*reinterpret_cast< ProblemParameters(*)>(_a[1]))); break;
        case 3: _t->quit(); break;
        case 4: _t->centerCamera(); break;
        case 5: _t->updateGL(); break;
        case 6: _t->makeCone(); break;
        case 7: _t->makeFlatCone(); break;
        case 8: _t->crush(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Controller::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Controller::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Controller,
      qt_meta_data_Controller, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Controller::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Controller::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Controller::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Controller))
        return static_cast<void*>(const_cast< Controller*>(this));
    return QObject::qt_metacast(_clname);
}

int Controller::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
