TEMPLATE = app
CONFIG += console c++16
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:/Libs/glm

SOURCES += \
    mesh.cpp \
    obj_import_export.cpp \
    mesh_motorcycle_graph.cpp \
    union_find.cpp \
    mesh_cutgraph_opt.cpp \
    primary_structure.cpp \
    primary_to_secondary.cpp \
    DNA.cpp \
    RNA.cpp \
    main.cpp
	
HEADERS += \
    union_find.h \
    primary_structure.h \
    RNA.h \
    mesh.h \
    DNA.h \
    uv.h

