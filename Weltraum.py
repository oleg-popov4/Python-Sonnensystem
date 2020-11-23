from scipy.integrate import solve_ivp, RK23
import numpy as np
import matplotlib.pyplot as plt
# import time


# Classen
class Planet(object):
    anzahl_planeten = 0

    def __init__(self, name):
        # Constructor
        self.pos = None
        self.geschw = None
        self.name = name  # String
        # 2 x n Matrix mit x und y Positionen vom Planeten
        # self.pos = pos
        # Anfangsbedingung fur die Berechnung
        # self.geschw = geschw
        Planet.anzahl_planeten = Planet.anzahl_planeten + 1


# End Planet

class Weltraum(Planet):
    # Konstanten
    Referenzzeit = 5.0 * 10.0 ** (6.0)                  # in  Sekunden
    Referenzzeit = Referenzzeit / (24.0 * 60.0 * 60.0)  # in  Tagen
    Referenzzeit_lange = 1495 * 10 ^ 8                  # Abstand zu Sonne

    # Funktionen
    def odefun1(self, t, vektor, komp_gros=4, dimension=2):
        epsilon = self.epsilon
        output = np.zeros(vektor.shape)
        anzah_planeten = int(min(vektor.size / komp_gros, len(epsilon)))
        pos_schrit = int(komp_gros / dimension)  # Brauche immer zwei Eintrage
        for iter1 in range(anzah_planeten):
            sum = 0.0
            pos1 = iter1 * komp_gros
            ri = vektor[pos1: pos1 + pos_schrit]
            dri = vektor[pos1 + pos_schrit: pos1 + komp_gros]
            # Setze dri - (Ableitung von ri) auf die neue Position, da dri nicht geandert wird
            output[pos1: pos1 + pos_schrit] = dri
            for iter2 in range(anzah_planeten):
                if not iter1 == iter2:
                    pos2 = iter2 * komp_gros
                    rj = vektor[pos2: pos2 + pos_schrit]
                    diferenz = ri - rj
                    sum = sum - epsilon[iter2] * diferenz / np.linalg.norm(diferenz) ** 3
                # Ende if
            # Ender for fuer iter2
            output[pos1 + pos_schrit: pos1 + komp_gros] = sum
        # Ende for fuer iter1
        return output

    # Ende def odefun1

    def weltraum_ploten(self):
        # Plote alle Planeten (Zeitschritte) in eine Funktion (Iterration)
        # Plot Einstelungen (Optionen) fuer Grafik
        plot_options1 = ['y--', 'w--', 'w--', 'w--', 'w--', 'w--', 'w--', 'w--', 'w--', 'w--', 'w--']
        plot_options2 = ['y', '#929591', '#ad8150', '#3f9b0b', 'r', '#e2ca76', '#ffffc2', '#ff028d']
        # setze background auf schwarz (Weltal)
        plt.style.use('dark_background')
        # !!! Включить интерактивный режим для анимации
        plt.ion()
        # Erzeuge eine Figure
        # plt.figure(str(time.time()))
        plt.figure()
        # Schleife uber alle Zeitschritte
        for itter in range(self.t_vektor.size):
            # !!! Очистить текущую фигуру
            plt.clf()
            # Отобразить график
            # Schleife ueber alle Planeten
            planet_namen = []
            for akt_planet, par1, par2 in zip(self.weltraum_planeten, plot_options1, plot_options2):
                planet_namen.append('Laufbahn der ' + str(akt_planet.name))
                planet_namen.append(str(akt_planet.name))
                # Gesamte laufbahn des Planeten
                plt.plot(akt_planet.pos[0, :], akt_planet.pos[1, :], par1, linewidth=0.1)
                # Aktuelle Position des Planeten
                plt.plot(akt_planet.pos[0, itter], akt_planet.pos[1, itter], marker='o', linestyle='None', color=par2)
            # End for
            plt.title('Time = '+ str(round(self.t_vektor[itter], 3)) + ', Erdtage = '
                      + str(round(self.t_vektor[itter]*Weltraum.Referenzzeit, 4)) + ' von '
                      + str(round(self.t_vektor[-1]*Weltraum.Referenzzeit, 4)))
            plt.legend(planet_namen, loc=0)
            plt.axis('equal')
            #plt.grid()
            # Установка отображаемых интервалов по осям
            # plt.xlim(0, maxSize)
            # plt.ylim(-1.1, 1.1)
            # !!! Следующие два вызова требуются для обновления графика
            plt.draw()
            plt.gcf().canvas.flush_events()
        # End for
        # Отключить интерактивный режим по завершению анимации
        plt.ioff()
        # Нужно, чтобы график не закрывался после завершения анимации
        # plt.show()

    # End Funktion

    def __init__(self, **kwargs):
        # Konstanten
        self.komp_gros = 4  # Feste grosse fur ein Planeten eintrag
        self.dimension = 2  # Planeten haben x und y Positionen
        # Default Startwerte
        # Default dictionary
        default_konstant = {'Tmax': 5, 'dt': 0.1, 'planeten_namen': ['Sonne', 'Erde', 'Jupiter'], 'komp_gros': 4,
                    'dimension': 2,
                    'anfangsbed_vektor': [np.array([0, 0]), np.array([0, 0]), np.array([1, 0]), np.array([0, 0.8]),
                                          np.array([5.2, 0]), np.array([0, 0.45])],
                    'epsilon' : [1, 3e-6, 9.6e-4]
                    }
        self.planeten_namen = kwargs.get('planeten_namen', default_konstant['planeten_namen'])
        self.anzah_planeten = len(self.planeten_namen)
        # Planeten position, geschwindigkeit des Planeten
        self.anfangsbed_vektor = kwargs.get('anfangsbed_vektor', default_konstant['anfangsbed_vektor'])
        self.epsilon = kwargs.get('epsilon', default_konstant['epsilon'])
        # Diskretisierung der Zeit
        self.Tmax = kwargs.get('Tmax', default_konstant['Tmax'])
        self.dt = kwargs.get('dt', default_konstant['dt'])
        self.N = self.Tmax / self.dt
        self.t_vektor = np.linspace(0, self.Tmax, int(self.N) + 1, endpoint=True)

        # Erstellung von Planeten
        self.weltraum_planeten = [Planet(iter1) for iter1 in self.planeten_namen]
        # Erstellung Startvektor
        self.vektor_start = np.array([], dtype=np.float64)
        for iter in self.anfangsbed_vektor:
            self.vektor_start = np.hstack((self.vektor_start, iter))
        # End for

        # End __init__

    def simulatn_start(self):
        print("Weltraum Simulation gestartet")
        # Berechnung der Loesung
        rhs = solve_ivp(self.odefun1, [0, self.Tmax], self.vektor_start, 'RK23', self.t_vektor, vectorized=True,
                        rtol=1e-10)
        # Eintrag der Loesungen
        for iter1, iter2 in zip(self.weltraum_planeten, range(self.anzah_planeten)):
            iter1.pos = rhs.y[iter2 * self.komp_gros: iter2 * self.komp_gros + self.dimension, :]
        # Loesungen Ploten
        self.weltraum_ploten()
        print("Weltraum Simulation beendet")
        return 1

# End Weltraum
