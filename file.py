import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def __init__(self, R=1, w=0.2, n=200):
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.X, self.Y, self.Z = self._compute_surface()
    
    def _compute_surface(self):
        U, V = self.U, self.V
        X = (self.R + V * np.cos(U / 2)) * np.cos(U)
        Y = (self.R + V * np.cos(U / 2)) * np.sin(U)
        Z = V * np.sin(U / 2)
        return X, Y, Z

    def surface_area(self):
        # Numerical integration using the magnitude of the cross product of partial derivatives
        du = 2 * np.pi / (self.n - 1)
        dv = self.w / (self.n - 1)
        
        Xu = np.gradient(self.X, du, axis=1)
        Yu = np.gradient(self.Y, du, axis=1)
        Zu = np.gradient(self.Z, du, axis=1)
        
        Xv = np.gradient(self.X, dv, axis=0)
        Yv = np.gradient(self.Y, dv, axis=0)
        Zv = np.gradient(self.Z, dv, axis=0)

        cross_prod_mag = np.sqrt((Yu * Zv - Zu * Yv) ** 2 +
                                 (Zu * Xv - Xu * Zv) ** 2 +
                                 (Xu * Yv - Yu * Xv) ** 2)
        
        area = np.sum(cross_prod_mag) * du * dv
        return area

    def edge_length(self):
        # Approximate the length along u at v = -w/2 and v = w/2
        def arc_length(Xc, Yc, Zc):
            dx = np.diff(Xc)
            dy = np.diff(Yc)
            dz = np.diff(Zc)
            return np.sum(np.sqrt(dx**2 + dy**2 + dz**2))

        edge1 = self.X[0, :], self.Y[0, :], self.Z[0, :]
        edge2 = self.X[-1, :], self.Y[-1, :], self.Z[-1, :]
        
        return arc_length(*edge1) + arc_length(*edge2)

    def plot(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.X, self.Y, self.Z, color='lightblue', edgecolor='k', alpha=0.7)
        ax.set_title("Mobius Strip")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.tight_layout()
        plt.show()

# Instantiate and compute properties
mobius = MobiusStrip(R=1, w=0.3, n=300)
surface_area = mobius.surface_area()
edge_length = mobius.edge_length()
mobius.plot()

(surface_area, edge_length)
