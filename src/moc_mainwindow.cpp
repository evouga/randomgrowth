/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon Apr 6 17:24:55 2015
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      28,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x0a,
      58,   52,   11,   11, 0x0a,
      97,   81,   11,   11, 0x0a,
     134,   11,   11,   11, 0x0a,
     148,   11,   11,   11, 0x0a,
     155,   11,   11,   11, 0x08,
     181,   11,   11,   11, 0x08,
     215,   11,   11,   11, 0x08,
     252,   11,   11,   11, 0x08,
     283,   11,   11,   11, 0x08,
     316,   11,   11,   11, 0x08,
     348,   11,   11,   11, 0x08,
     385,  380,   11,   11, 0x08,
     425,  380,   11,   11, 0x08,
     466,  380,   11,   11, 0x08,
     503,  380,   11,   11, 0x08,
     538,  380,   11,   11, 0x08,
     578,  380,   11,   11, 0x08,
     619,  380,   11,   11, 0x08,
     657,  380,   11,   11, 0x08,
     697,  380,   11,   11, 0x08,
     730,  380,   11,   11, 0x08,
     764,  380,   11,   11, 0x08,
     801,  380,   11,   11, 0x08,
     844,   11,   11,   11, 0x08,
     879,   11,   11,   11, 0x08,
     910,   11,   11,   11, 0x08,
     935,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0params\0setParameters(ProblemParameters)\0"
    "error\0showError(std::string)\0"
    "centroid,radius\0centerCamera(Eigen::Vector3d,double)\0"
    "repaintMesh()\0tick()\0on_actionExit_triggered()\0"
    "on_actionReset_Camera_triggered()\0"
    "on_actionTake_Screenshot_triggered()\0"
    "on_wireframeCheckBox_clicked()\0"
    "on_smoothShadeCheckBox_clicked()\0"
    "on_actionExport_OBJ_triggered()\0"
    "on_actionImport_OBJ_triggered()\0arg1\0"
    "on_poissonRatioEdit_textEdited(QString)\0"
    "on_youngsModulusEdit_textEdited(QString)\0"
    "on_thicknessEdit_textEdited(QString)\0"
    "on_densityEdit_textEdited(QString)\0"
    "on_dampingCoeffEdit_textEdited(QString)\0"
    "on_eulerTimestepEdit_textEdited(QString)\0"
    "on_eulerItersEdit_textEdited(QString)\0"
    "on_growthAmountEdit_textEdited(QString)\0"
    "on_scaleEdit_textEdited(QString)\0"
    "on_outputEdit_textEdited(QString)\0"
    "on_maxStrainEdit_textEdited(QString)\0"
    "on_baseProbabilityEdit_textEdited(QString)\0"
    "on_actionImport_Metric_triggered()\0"
    "on_actionMake_Cone_triggered()\0"
    "on_crushButton_clicked()\0"
    "on_actionMake_Flat_Cone_triggered()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->setParameters((*reinterpret_cast< ProblemParameters(*)>(_a[1]))); break;
        case 1: _t->showError((*reinterpret_cast< std::string(*)>(_a[1]))); break;
        case 2: _t->centerCamera((*reinterpret_cast< Eigen::Vector3d(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 3: _t->repaintMesh(); break;
        case 4: _t->tick(); break;
        case 5: _t->on_actionExit_triggered(); break;
        case 6: _t->on_actionReset_Camera_triggered(); break;
        case 7: _t->on_actionTake_Screenshot_triggered(); break;
        case 8: _t->on_wireframeCheckBox_clicked(); break;
        case 9: _t->on_smoothShadeCheckBox_clicked(); break;
        case 10: _t->on_actionExport_OBJ_triggered(); break;
        case 11: _t->on_actionImport_OBJ_triggered(); break;
        case 12: _t->on_poissonRatioEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 13: _t->on_youngsModulusEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: _t->on_thicknessEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: _t->on_densityEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 16: _t->on_dampingCoeffEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 17: _t->on_eulerTimestepEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 18: _t->on_eulerItersEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 19: _t->on_growthAmountEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 20: _t->on_scaleEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 21: _t->on_outputEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 22: _t->on_maxStrainEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 23: _t->on_baseProbabilityEdit_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 24: _t->on_actionImport_Metric_triggered(); break;
        case 25: _t->on_actionMake_Cone_triggered(); break;
        case 26: _t->on_crushButton_clicked(); break;
        case 27: _t->on_actionMake_Flat_Cone_triggered(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 28)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 28;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
