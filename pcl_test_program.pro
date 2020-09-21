QT += gui

CONFIG += c++11
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        alglib/alglibinternal.cpp \
        alglib/alglibmisc.cpp \
        alglib/ap.cpp \
        alglib/dataanalysis.cpp \
        alglib/diffequations.cpp \
        alglib/fasttransforms.cpp \
        alglib/integration.cpp \
        alglib/interpolation.cpp \
        alglib/linalg.cpp \
        alglib/optimization.cpp \
        alglib/solvers.cpp \
        alglib/specialfunctions.cpp \
        alglib/statistics.cpp \
        main.cpp \
        mathstatistics.cpp

HEADERS += \
        alglib/alglibinternal.h \
        alglib/alglibmisc.h \
        alglib/ap.h \
        alglib/dataanalysis.h \
        alglib/diffequations.h \
        alglib/fasttransforms.h \
        alglib/integration.h \
        alglib/interpolation.h \
        alglib/linalg.h \
        alglib/optimization.h \
        alglib/solvers.h \
        alglib/specialfunctions.h \
        alglib/statistics.h \
        alglib/stdafx.h \
        mathstatistics.h \
        stl_reader.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

INCLUDEPATH = $$PWD/alglib

include(pcl_library_path.pri)
include(opencascade_library_path.pri)
