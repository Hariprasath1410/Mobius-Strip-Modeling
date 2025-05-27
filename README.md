# Mobius-Strip-Modeling
Code Structure & Explanation
1. Class Design:

MobiusStrip: Encapsulates the strip parameters, mesh generation, and geometry calculations.

_compute_surface(): Implements parametric equations using numpy.meshgrid.

surface_area(): Approximates the area using a vector cross product of partial derivatives, mimicking surface integral via Riemann sums.

edge_length(): Computes arc length of two boundary edges using Euclidean distances.

plot(): Visualizes the 3D Mobius strip.

Outputs
Surface Area ≈ 1.90 units²

Edge Length ≈ 12.60 units

Challenges
Capturing the twist correctly in discretized grid.

Ensuring smooth derivative estimation for accurate area approximation.

Visualizing the non-orientable nature clearly with only 2D projections.
