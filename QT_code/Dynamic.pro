#-------------------------------------------------
#
# Project created by QtCreator 2018-08-13T15:00:33
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Dynamic
TEMPLATE = app


SOURCES += main.cpp\
        dynamicwidget.cpp \
    c_algorithm/algo.c \
    c_algorithm/dynIdenf.c \
    c_algorithm/dynIdenf_emxAPI.c \
    c_algorithm/dynIdenf_emxutil.c \
    c_algorithm/loadfile.c \
    c_algorithm/write_results.c \
    c_algorithm/write_error.c

HEADERS  += dynamicwidget.h \
    c_algorithm/dynIdenf.h \
    c_algorithm/dynIdenf_emxAPI.h \
    c_algorithm/dynIdenf_emxutil.h \
    c_algorithm/dynIdenf_types.h \
    c_algorithm/rtwtypes.h

FORMS    += dynamicwidget.ui

RESOURCES += \
    image.qrc

DISTFILES +=

CONFIG += C++11
QMAKE_CFLAGS += -std=c99
