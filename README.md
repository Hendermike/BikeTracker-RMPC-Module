# Bike Tracker's RMPC module

This module is one of the bulding blocks of the Bike Tracker project. Programmed on Java to enable its direct integration with Android Studio projects.<br />

Implementation of a robust model-based predictive control module for cyclists acceleration and fuzzy acceleration-error intervals calculation, required to enable platoon driving behavior (constant distance between platoon members and same speed as the leader).
Required accleration predictions calculation (for 5 future samples) is formulated as a Adaptive Cruise Control problem, solved using a genetic algorithm strategy on a ad-hoc cuadratic functional based optimization problem.
Robust intervals for the acceleration-error are estimated using a Takagi-Sugeno model built upon human acceleration-error historical data.

