import numpy as np
from Weltraum import Weltraum

if __name__ == '__main__':
    Referenzzeit = 5.0 * 10.0 ** (6.0) / (24.0 * 60.0 * 60.0)  # in  Tagen
    anzahl_Planeten = 5
    Tmax = 1*366/Referenzzeit
    dt = 0.5/Referenzzeit
    planeten_namen = ['Sonne', 'Merkur', 'Venus', 'Erde', 'Mars', 'Jupiter', 'Saturn']
    #                  | position  x                   y      |geschwind x  y |                     Planet
    anfangsbed_vektor= [np.array([0,                  0]),     np.array([0, 0]),                    # Sonne
                        np.array([0.387959866220736,  0]),     np.array([0, 1.58651072183211]),     # Merkur
                        np.array([0.722408026755853,  0]),     np.array([0, 1.17313356162501]),     # Venus
                        np.array([1,                  0]),     np.array([0, 1.005]),                # Erde
                        np.array([1.52508361204013,   0]),     np.array([0, 0.806319955120330]),    # Mars
                        np.array([5.20401337792642,   0]),     np.array([0, 0.437496411045763]),    # Jupiter
                        np.array([9.53846153846154, 0]),       np.array([0, 0.324269927942035]),    # Saturn
                        ]
    #         'Sonne','Merkur',         'Venus',        'Erde',         'Mars',        'Jupiter'    'Saturn'
    epsilon = [1,     1.66013*10**-7, 2.45172*10**-6, 3.00362*10**-6, 3.22722*10**-7, 0.000955039, 0.000285807684570509
               ]

    simulation1 = Weltraum(Tmax=Tmax, dt=dt, planeten_namen=planeten_namen[0:anzahl_Planeten], anfangsbed_vektor=anfangsbed_vektor[0:anzahl_Planeten*4],
                           epsilon=epsilon[0:anzahl_Planeten])
    simulation1.simulatn_start()
