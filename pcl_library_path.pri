#win32:CONFIG(release, debug|release): LIBS += -L'C:/Program Files/PCL 1.11.1/lib/' -lpcl_common \
LIBS += -L'C:/Program Files/PCL 1.11.1/lib/' -lpcl_common \
    -lpcl_features \
    -lpcl_filters \
    -lpcl_io \
    -lpcl_io_ply \
    -lpcl_kdtree \
    -lpcl_keypoints \
    -lpcl_ml \
    -lpcl_octree \
    -lpcl_outofcore \
    -lpcl_people \
    -lpcl_recognition \
    -lpcl_registration \
    -lpcl_sample_consensus \
    -lpcl_search \
    -lpcl_segmentation \
    -lpcl_stereo \
    -lpcl_surface \
    -lpcl_tracking \
    -lpcl_visualization

#else:win32:CONFIG(debug, debug|release): LIBS += -L'C:/Program Files/PCL 1.11.1/lib/' -lpcl_commond \
#    -lpcl_featuresd \
#    -lpcl_filtersd \
#    -lpcl_iod \
#    -lpcl_io_plyd \
#    -lpcl_kdtreed \
#    -lpcl_keypointsd \
#    -lpcl_mld \
#    -lpcl_octreed \
#    -lpcl_outofcored \
#    -lpcl_peopled \
#    -lpcl_recognitiond \
#    -lpcl_registrationd \
#    -lpcl_sample_consensusd \
#    -lpcl_searchd \
#    -lpcl_segmentationd \
#    -lpcl_stereod \
#    -lpcl_surfaced \
#    -lpcl_trackingd \
#    -lpcl_visualizationd

INCLUDEPATH += 'C:\Program Files\PCL 1.11.1\include\pcl-1.11'
DEPENDPATH += 'C:\Program Files\PCL 1.11.1\include\pcl-1.11'

include(eigen_library_path.pri)
include(boots_library_path.pri)
include(flann_library_path.pri)
include(vtk_library_path.pri)
