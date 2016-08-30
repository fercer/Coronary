/****************************************************************************
** Meta object code from reading C++ file 'coronarygui.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "coronarygui.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'coronarygui.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_coronaryGUI[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      41,   28,   12,   12, 0x08,
      72,   12,   12,   12, 0x08,
      83,   12,   12,   12, 0x08,
     101,   12,   12,   12, 0x08,
     133,   12,   12,   12, 0x08,
     160,   12,   12,   12, 0x08,
     193,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_coronaryGUI[] = {
    "coronaryGUI\0\0defineColors()\0,viewport_ID\0"
    "showImage(IMGVTK::IMG_IDX,int)\0"
    "loadBase()\0loadGroundtruth()\0"
    "on_action_Open_file_triggered()\0"
    "on_action_Quit_triggered()\0"
    "on_actionSkeletonize_triggered()\0"
    "on_actionOpen_file_as_ground_truth_triggered()\0"
};

void coronaryGUI::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        coronaryGUI *_t = static_cast<coronaryGUI *>(_o);
        switch (_id) {
        case 0: _t->defineColors(); break;
        case 1: _t->showImage((*reinterpret_cast< IMGVTK::IMG_IDX(*)>(_a[1])),(*reinterpret_cast< const int(*)>(_a[2]))); break;
        case 2: _t->loadBase(); break;
        case 3: _t->loadGroundtruth(); break;
        case 4: _t->on_action_Open_file_triggered(); break;
        case 5: _t->on_action_Quit_triggered(); break;
        case 6: _t->on_actionSkeletonize_triggered(); break;
        case 7: _t->on_actionOpen_file_as_ground_truth_triggered(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData coronaryGUI::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject coronaryGUI::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_coronaryGUI,
      qt_meta_data_coronaryGUI, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &coronaryGUI::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *coronaryGUI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *coronaryGUI::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_coronaryGUI))
        return static_cast<void*>(const_cast< coronaryGUI*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int coronaryGUI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
