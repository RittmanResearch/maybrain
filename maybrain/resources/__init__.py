"""
Initialisation file for resources package

"""
import pkg_resources

PROPERTIES_ANATLABEL_500 = pkg_resources.resource_filename('maybrain', 'resources/properties_anatlabels_500.txt')
PROPERTIES_HEMISPHERES_500 = pkg_resources.resource_filename('maybrain', 'resources/properties_hemispheres_500.txt')
PROPERTIES_LOBES_500 = pkg_resources.resource_filename('maybrain', 'resources/properties_lobes_500.txt')

DUMMY_ADJ_FILE_500 = pkg_resources.resource_filename('maybrain', 'resources/adj_file_500.txt')

MNI_SPACE_COORDINATES_500 = pkg_resources.resource_filename('maybrain', 'resources/parcel_500.txt')
