#win32:CONFIG(release, debug|release): LIBS += -L'C:/Program Files/PCL 1.11.1/3rdParty/FLANN/lib/' -lflann \
LIBS += -L'C:/Program Files/PCL 1.11.1/3rdParty/FLANN/lib/' -lflann \
    -lflann_cpp \
    -lflann_cpp_s \
    -lflann_s
#else:win32:CONFIG(debug, debug|release): LIBS += -L'C:/Program Files/PCL 1.11.1/3rdParty/FLANN/lib/' -lflann-gd \
#    -lflann_cpp-gd \
#    -lflann_cpp_s-gd \
#    -lflann_s-gd \

INCLUDEPATH += 'C:\Program Files\PCL 1.11.1\3rdParty\FLANN\include'
DEPENDPATH += 'C:\Program Files\PCL 1.11.1\3rdParty\FLANN\include'
