"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


"""default values for each data type"""

default_values_HF = {
    'Downsampling rate': '30',
    'Butterworth filter degree': '4',
    'Cut-off frequency': '0.005',
    'Force threshold, pN': '5',
    'Z-score force': '3',
    'Z-score distance': '3',
    'Step d': '10',
    'Moving median window size': '800',
    'STD difference threshold': '0.05',
    'Data frequency, Hz': '1000'
}

default_values_LF = {
    'Downsampling rate': '1',
    'Butterworth filter degree': '2',
    'Cut-off frequency': '0.5',
    'Force threshold, pN': '5',
    'Z-score force': '3',
    'Z-score distance': '3',
    'Step d': '3',
    'Moving median window size': '20',
    'STD difference threshold': '0.05',
    'Data frequency, Hz': '1000'
}

default_values_CSV = {
    'Downsampling rate': '1',
    'Butterworth filter degree': '1',
    'Cut-off frequency': '0.01',
    'Force threshold, pN': '5',
    'Z-score force': '2.5',
    'Z-score distance': '3',
    'Step d': '10',
    'Moving median window size': '120',
    'STD difference threshold': '0.05',
    'Data frequency, Hz': '20'
}

default_values_FIT = {
    'Persistance-Length ds, nm': '40',
    'Persistance-Length ds, upper bound, nm': '80',
    'Persistance-Length ds, lower bound, nm': '12',
    'Persistance-Length ss, nm': '1',
    'Contour-Length ds, nm': '1260',
    'Contour-Length ss, nm': '1',
    'Contour-Length ss, upper bound, nm': '1267',
    'Stiffness ds, pN': '720',
    'Stiffness ds, upper bound, pN': '750',
    'Stiffness ds, lower bound, pN': '700',
    'Stiffness ss, pN': '1000',
    'Stiffness ss, upper bound, pN': '1100',
    'Stiffness ss, lower bound, pN': '900',
    'Force offset, pN': '0',
    'Force offset, upper bound, pN': '1',
    'Force offset, lower bound, pN': '-1',
    'Distance offset, nm': '0',
    'Distance offset, upper bound, nm': '500',
    'Distance offset, lower bound, nm': '-500'
}

default_values_constantF = {
    'x min': '0',
    'x max': '180',
    'y min': '1290',
    'y max': '1320',
    'Number gauss': '3',
    'Mean': '1295,1310,1317',
    'STD': '1,1,1',
    'Amplitude': '10000,5000,10000'
}
