/****************************************************************************
** Meta object code from reading C++ file 'filterpars.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "filterpars.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'filterpars.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_filterpars[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      20,   12,   11,   11, 0x08,
      45,   12,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_filterpars[] = {
    "filterpars\0\0checked\0on_chkTune_toggled(bool)\0"
    "on_chkLFix_toggled(bool)\0"
};

void filterpars::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        filterpars *_t = static_cast<filterpars *>(_o);
        switch (_id) {
        case 0: _t->on_chkTune_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->on_chkLFix_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData filterpars::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject filterpars::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_filterpars,
      qt_meta_data_filterpars, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &filterpars::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *filterpars::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *filterpars::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_filterpars))
        return static_cast<void*>(const_cast< filterpars*>(this));
    return QDialog::qt_metacast(_clname);
}

int filterpars::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
